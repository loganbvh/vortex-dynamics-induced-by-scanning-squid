import json
import logging
import os
import sys
from datetime import datetime
from typing import Any, Dict, Optional, Sequence, Tuple, Union

import h5py
import numpy as np
import superscreen as sc
import tdgl
from tdgl.geometry import box, circle
from tdgl.solution.data import get_data_range

sys.path.append("..")
import squids

squid_funcs = {
    "ibm-small": squids.ibm.small.make_squid,
    "ibm-medium": squids.ibm.medium.make_squid,
    "ibm-large": squids.ibm.large.make_squid,
    "ibm-xlarge": squids.ibm.xlarge.make_squid,
    "huber": squids.huber.make_squid,
    "hypres-small": squids.hypres.small.make_squid,
}


logger = logging.getLogger(os.path.basename(__file__).replace(".py", ""))
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler(sys.stdout))


def make_film(
    radius: float,
    xi: float,
    lambda_: float,
    d: float,
    min_points: int,
    smooth: int,
    max_edge_length: Optional[float] = None,
    shape: str = "circle",
    points: int = 501,
) -> tdgl.Device:
    """Make the tdgl.Device representing the Nb film."""
    layer = tdgl.Layer(coherence_length=xi, london_lambda=lambda_, thickness=d, z0=0)
    assert shape in ("box", "circle"), shape
    if shape == "circle":
        film = tdgl.Polygon("film", points=circle(radius, points=points))
    else:
        film = tdgl.Polygon("film", points=box(2 * radius, points=points)).buffer(0)
    device = tdgl.Device(
        shape,
        layer=layer,
        film=film,
        length_units="um",
    )
    device.make_mesh(
        min_points=min_points, max_edge_length=max_edge_length, smooth=smooth
    )
    print(device.mesh_stats_dict())
    return device


def make_film_with_slot(
    radius: float,
    xi: float,
    lambda_: float,
    d: float,
    min_points: int,
    smooth: int,
    slot_size: Tuple[float, float],
    slot_top_center: Tuple[float, float],
    slot_radius: float = 0,
    max_edge_length: Optional[float] = None,
    shape: str = "circle",
    points: int = 501,
) -> tdgl.Device:
    """Make the tdgl.Device representing the Nb film."""
    layer = tdgl.Layer(coherence_length=xi, london_lambda=lambda_, thickness=d, z0=0)
    assert shape in ("box", "circle"), shape
    if shape == "circle":
        film = tdgl.Polygon("film", points=circle(radius, points=points))
    else:
        film = tdgl.Polygon("film", points=box(2 * radius, points=points)).buffer(0)
    width, height = np.array(slot_size) - 2 * slot_radius
    dx, dy = np.array(slot_top_center)
    slot = (
        tdgl.Polygon("slot", points=box(width, height, points=201))
        .translate(dy=-(height + slot_radius) / 2)
        .translate(dx=dx, dy=dy)
    )
    if slot_radius:
        slot = slot.buffer(slot_radius, join_style="round").resample(201)
    if film.contains_points(slot.points).all():
        holes = [slot]
    else:
        film = film.difference(slot).resample(501)
        holes = None
    device = tdgl.Device(
        "box",
        layer=layer,
        film=film,
        holes=holes,
        length_units="um",
    )
    device.make_mesh(
        min_points=min_points, max_edge_length=max_edge_length, smooth=smooth
    )
    print(device.mesh_stats_dict())
    return device


def make_squid(
    squid_type: str, min_points: int, smooth: int, angle: float
) -> sc.Device:
    """Make the superscreen.Device representing the SQUID."""
    squid = squid_funcs[squid_type]().rotate(angle)
    squid.make_mesh(min_points=min_points, smooth=smooth)
    return squid


def get_base_squid_solution(squid: sc.Device, iterations: int) -> sc.Solution:
    """Generate the supercreen.Solution for 1 mA in the SQUID field coil."""
    return sc.solve(
        squid,
        circulating_currents=dict(fc_center=1.0),
        current_units="mA",
        field_units="mT",
        iterations=iterations,
    )[-1]


def applied_potential(
    x: Sequence[float],
    y: Sequence[float],
    z: Union[float, Sequence[float]],
    *,
    path_to_solution: str,
    I_fc: float,
    r0: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    current_units: str = "mA",
    field_units: str = "mT",
    length_units: str = "um",
) -> np.ndarray:
    """Calculates the vector potential applied to the film for a
    given current in the SQUID field coil.
    """
    solution = sc.Solution.from_file(path_to_solution, compute_matrices=True)
    r0 = np.atleast_2d(r0)
    if len(z) == 1:
        z = z[0] * np.ones_like(x)
    positions = np.array([x.squeeze(), y.squeeze(), z.squeeze()]).T - r0
    I_circ = solution.circulating_currents["fc_center"] * tdgl.ureg(
        solution.current_units
    )
    I_fc = I_fc * tdgl.ureg(current_units)
    scale = (I_fc / I_circ).to_base_units().magnitude
    A = solution.vector_potential_at_position(
        positions,
        units=f"{field_units} * {length_units}",
        with_units=False,
    )
    return scale * A


def applied_field_from_film(
    x: Sequence[float],
    y: Sequence[float],
    z: Union[float, Sequence[float]],
    *,
    tdgl_solution: tdgl.Solution,
    field_units: str = "mT",
    with_units: bool = False,
):
    """Calculates the magnetic field due to currents flowing in the film."""
    if isinstance(z, (int, float)):
        z = [z]
    if len(z) == 1:
        z = z[0] * np.ones_like(x)
    positions = np.array([x.squeeze(), y.squeeze(), z.squeeze()]).T
    return tdgl_solution.field_at_position(
        positions,
        vector=False,
        units=field_units,
        with_units=with_units,
        return_sum=True,
    )


def simulate_dynamics_squid(
    *,
    directory: str,
    device: tdgl.Device,
    squid: sc.Device,
    squid_position: Tuple[float, float, float],
    squid_iterations: int,
    field_units: str,
    current_units: str,
    I_fc: Sequence[float],
    seed_solutions: bool,
    solve_time: Tuple[float, float],
    dt_init: float,
    screening: bool,
    save_every: int,
    metadata: Dict[str, Any],
) -> None:
    """Simulates a SQUID suceptibility measurement for a given set
    of field coil currents.
    """

    start_time = datetime.now()
    path = os.path.abspath(directory)
    os.makedirs(path, exist_ok=True)

    squid_solution = get_base_squid_solution(squid, squid_iterations)
    squid_solution_path = os.path.join(path, "squid_solution")
    squid_solution.to_file(squid_solution_path)

    field_units = field_units
    current_units = current_units
    length_units = "um"

    all_flux = []
    all_fluxoids = []

    def calculate_pl_flux(tdgl_solution: tdgl.Solution) -> float:
        """Calculates the flux through the pickup loop due to currents flowing in the film."""
        applied_field = sc.Parameter(
            applied_field_from_film,
            tdgl_solution=tdgl_solution,
            field_units="mT",
        )
        solution = sc.solve(
            squid,
            applied_field=applied_field,
            iterations=squid_iterations,
            field_units="mT",
            current_units="mA",
        )[-1]
        fluxoid = solution.hole_fluxoid("pl_center")
        return sum(fluxoid).to("Phi_0").magnitude

    prev_solution = None

    with h5py.File(
        os.path.join(path, "steady-state.h5"),
        "x",
        libver="latest",
    ) as f:
        device.to_hdf5(f.create_group("solution/device"))

    for i, current in enumerate(I_fc):
        options = tdgl.SolverOptions(
            solve_time=(solve_time[1] if i and seed_solutions else solve_time[0]),
            dt_init=dt_init,
            output_file=os.path.join(path, f"output-{i}.h5"),
            save_every=save_every,
            include_screening=screening,
            field_units=field_units,
        )

        A_applied = tdgl.Parameter(
            applied_potential,
            path_to_solution=squid_solution_path,
            r0=squid_position,
            I_fc=current,
            field_units=field_units,
            current_units=current_units,
            length_units=length_units,
        )

        tdgl_solution = tdgl.solve(
            device,
            options,
            applied_vector_potential=A_applied,
            seed_solution=prev_solution,
        )
        tdgl_solution.to_hdf5()

        if seed_solutions:
            prev_solution = tdgl_solution

        flux = calculate_pl_flux(tdgl_solution)
        all_flux.append(flux)

        total_phase = tdgl_solution.boundary_phases(delta=True)["film"].phases[-1]
        all_fluxoids.append(total_phase / (2 * np.pi))

        with h5py.File(tdgl_solution.path, "r+", libver="latest") as f:
            for key, val in metadata.items():
                if val is not None:
                    f.attrs[key] = val
            f.attrs["pl_fluxoid_in_Phi_0"] = flux

            _, i_end = get_data_range(f)

            with h5py.File(os.path.join(path, "steady-state.h5"), "r+") as out:
                data_grp = out.require_group("data")
                f["data"].copy(str(i_end), data_grp, name=str(i))
                for key, val in metadata.items():
                    if val is not None:
                        out.attrs[key] = val
                data_grp[str(i)].attrs["pl_fluxoid_in_Phi_0"] = flux

    end_time = datetime.now()

    json_data = {}
    json_data["args"] = metadata.copy()
    json_data["timing"] = {
        "start_time": start_time.isoformat(),
        "end_time": end_time.isoformat(),
        "total_seconds": (end_time - start_time).total_seconds(),
    }
    json_data["I_fc"] = I_fc.tolist()
    json_data["flux"] = all_flux
    json_data["film_fluxoid"] = all_fluxoids

    with open(os.path.join(path, "results.json"), "w") as f:
        json.dump(json_data, f, indent=4, sort_keys=True)


def simulate_dynamics_loop(
    *,
    directory: str,
    device: tdgl.Device,
    fc_radius: float,
    fc_center: Tuple[float, float, float],
    pl_radius: float,
    pl_center: Tuple[float, float, float],
    field_units: str,
    current_units: str,
    I_fc: Sequence[float],
    seed_solutions: bool,
    solve_time: Tuple[float, float],
    dt_init: float,
    screening: bool,
    save_every: int,
    metadata: Dict[str, Any],
) -> None:
    """Simulates a SQUID suceptibility measurement for a given set
    of field coil currents.
    """

    start_time = datetime.now()
    path = os.path.abspath(directory)
    os.makedirs(path, exist_ok=True)

    field_units = field_units
    current_units = current_units
    length_units = "um"

    all_flux = []
    all_fluxoids = []

    pl = tdgl.Polygon(points=circle(pl_radius, center=pl_center[:2], points=201))
    pl_mesh = pl.make_mesh(min_points=500, smooth=100)

    def calculate_pl_flux(tdgl_solution):
        """Calculates the flux through the pickup loop due to currents flowing in the film."""
        fields = applied_field_from_film(
            pl_mesh.x,
            pl_mesh.y,
            pl_center[-1],
            tdgl_solution=tdgl_solution,
            with_units=True,
        )
        areas = pl_mesh.areas * tdgl.ureg(length_units) ** 2
        return np.sum(fields * areas).to("Phi_0").magnitude

    prev_solution = None

    with h5py.File(
        os.path.join(path, "steady-state.h5"),
        "x",
        libver="latest",
    ) as f:
        device.to_hdf5(f.create_group("solution/device"))

    for i, current in enumerate(I_fc):
        options = tdgl.SolverOptions(
            solve_time=(solve_time[1] if i and seed_solutions else solve_time[0]),
            dt_init=dt_init,
            output_file=os.path.join(path, f"output-{i}.h5"),
            save_every=save_every,
            include_screening=screening,
            field_units=field_units,
        )

        A_applied = tdgl.sources.CurrentLoop(
            current=current,
            radius=fc_radius,
            center=fc_center,
            current_units=current_units,
            field_units=field_units,
            length_units=length_units,
        )

        tdgl_solution = tdgl.solve(
            device,
            options,
            applied_vector_potential=A_applied,
            seed_solution=prev_solution,
        )
        tdgl_solution.to_hdf5()

        if seed_solutions:
            prev_solution = tdgl_solution

        flux = calculate_pl_flux(tdgl_solution)
        all_flux.append(flux)

        total_phase = tdgl_solution.boundary_phases(delta=True)["film"].phases[-1]
        all_fluxoids.append(total_phase / (2 * np.pi))

        with h5py.File(tdgl_solution.path, "r+", libver="latest") as f:
            for key, val in metadata.items():
                if val is not None:
                    f.attrs[key] = val
            f.attrs["pl_fluxoid_in_Phi_0"] = flux

            _, i_end = get_data_range(f)

            with h5py.File(os.path.join(path, "steady-state.h5"), "r+") as out:
                data_grp = out.require_group("data")
                f["data"].copy(str(i_end), data_grp, name=str(i))
                for key, val in metadata.items():
                    if val is not None:
                        out.attrs[key] = val
                data_grp[str(i)].attrs["pl_fluxoid_in_Phi_0"] = flux

    end_time = datetime.now()

    json_data = {}
    json_data["args"] = metadata.copy()
    json_data["timing"] = {
        "start_time": start_time.isoformat(),
        "end_time": end_time.isoformat(),
        "total_seconds": (end_time - start_time).total_seconds(),
    }
    json_data["I_fc"] = I_fc.tolist()
    json_data["flux"] = all_flux
    json_data["film_fluxoid"] = all_fluxoids

    with open(os.path.join(path, "results.json"), "w") as f:
        json.dump(json_data, f, indent=4, sort_keys=True)


def main():
    import argparse

    parser = argparse.ArgumentParser()

    sample_grp = parser.add_argument_group("sample")
    squid_grp = parser.add_argument_group("squid")
    tdgl_grp = parser.add_argument_group("tdgl")

    squid_grp.add_argument(
        "--squid-type",
        type=str,
        default="hypres-small",
        choices=list(squid_funcs) + ["loop"],
    )
    squid_grp.add_argument(
        "--squid-points",
        type=int,
        default=4000,
        help="Minimum number of points in the SQUID mesh.",
    )
    squid_grp.add_argument(
        "--squid-smooth",
        type=int,
        default=50,
        help="Number of mesh smoohting steps for the SQUID mesh.",
    )
    squid_grp.add_argument(
        "--squid-angle",
        type=float,
        default=0,
        help="Angle by which to rotate the SQUID device, in degrees.",
    )
    squid_grp.add_argument(
        "--squid-position",
        type=float,
        nargs=3,
        default=0,
        help="SQUID (x, y, z) position in microns.",
    )
    squid_grp.add_argument(
        "--squid-iterations",
        type=int,
        default=5,
        help="Number of superscreen solve iterations.",
    )
    squid_grp.add_argument(
        "--fc-radius",
        type=float,
        help="1D loop field coil radius",
        default=2.4,
    )
    squid_grp.add_argument(
        "--fc-center",
        type=float,
        nargs=3,
        default=(0, 0, 1),
        help="1D loop field coil center position.",
    )
    squid_grp.add_argument(
        "--pl-radius",
        type=float,
        help="1D pickup loop radius",
        default=0.8,
    )
    squid_grp.add_argument(
        "--pl-center",
        type=float,
        nargs=3,
        default=None,
        help="1D pickup loop center position.",
    )

    sample_grp.add_argument(
        "--film-radius", type=float, default=15, help="Film radius in microns."
    )
    sample_grp.add_argument(
        "--film-shape",
        type=str,
        choices=("circle", "box"),
        default="circle",
        help="Shape of the film",
    )
    sample_grp.add_argument(
        "--slot-size",
        type=float,
        nargs=2,
        default=None,
        help="Width and height of the slot, if any.",
    )
    sample_grp.add_argument(
        "--slot-top-center-x",
        type=float,
        nargs="*",
        default=None,
        help="The x position of the top center of the slot.",
    )
    sample_grp.add_argument(
        "--slot-top-center-y",
        type=float,
        nargs="*",
        default=None,
        help="The y position of the top center of the slot.",
    )
    sample_grp.add_argument(
        "--slot-radius",
        type=float,
        default=0,
        help="Slot radius of curvature in microns.",
    )
    sample_grp.add_argument(
        "--film-points",
        type=int,
        default=4000,
    )
    sample_grp.add_argument(
        "--film-smooth",
        type=int,
        default=100,
    )
    sample_grp.add_argument(
        "--max-edge-length",
        type=float,
        default=0.5,
        help="Maximum edge length in the film mesh.",
    )
    sample_grp.add_argument(
        "--d", default=0.1, type=float, help="Film thickness in microns."
    )
    sample_grp.add_argument(
        "--lam",
        default=2,
        type=float,
        help="London penetration depth in microns.",
    )
    sample_grp.add_argument(
        "--xi",
        default=1,
        type=float,
        help="Coherence length in microns.",
    )
    sample_grp.add_argument(
        "--gamma", default=10, type=float, help="TDGL gamma parameter."
    )
    sample_grp.add_argument(
        "--screening", action="store_true", help="Include screening."
    )

    tdgl_grp.add_argument("--directory", type=str, help="Output directory.")
    tdgl_grp.add_argument(
        "--I_fc",
        nargs="*",
        type=float,
        help="Peak field coil current in mA: start, stop, num_steps.",
    )
    tdgl_grp.add_argument(
        "--index",
        type=int,
        help="Index for peak current in I_fc array.",
    )
    tdgl_grp.add_argument(
        "--cycles",
        default=3,
        type=float,
        help="Number of AC field cycles.",
    )
    tdgl_grp.add_argument(
        "--points-per-cycle",
        type=float,
        default=10,
        help="Number of current points per AC cycle.",
    )
    tdgl_grp.add_argument(
        "--seed-solutions",
        action="store_true",
        help="Seed each simulation with the previous solution.",
    )
    tdgl_grp.add_argument(
        "--solve-time",
        type=float,
        nargs=2,
        help="Solve time in units of GL tau.",
    )
    tdgl_grp.add_argument(
        "--dt-init",
        type=float,
        default=1e-6,
        help="Initial time step.",
    )
    tdgl_grp.add_argument(
        "--save-every", default=100, type=int, help="Save interval in steps."
    )
    tdgl_grp.add_argument(
        "--field-units",
        type=str,
        default="mT",
    )
    tdgl_grp.add_argument(
        "--current-units",
        type=str,
        default="mA",
    )

    args = parser.parse_args()
    args_as_dict = vars(args)
    for k, v in args_as_dict.items():
        print(f"{k}: {v}")

    if args.slot_size is None:
        print("No slot defined. Sweeping peak field coil current.")
        if len(args.I_fc) != 3:
            raise ValueError("I_fc must be a 3-tuple, (start, stop, num).")
        start, stop, num = args.I_fc
        I_fc_pk = np.linspace(start, stop, int(num))[args.index]
    elif len(args.I_fc) > 1:
        print("Slot exists. Sweeping peak field coil current.")
        start, stop, num = args.I_fc[:3]
        I_fc_pk = np.linspace(start, stop, int(num))[args.index]
        slot_top_center = (args.slot_top_center_x[0], args.slot_top_center_y[0])
    else:
        print("Slot exists. Sweeping slot position at fixed peak field coil current.")
        if len(args.I_fc) > 1:
            raise ValueError("I_fc must be a 1-tuple, (I_fc_pk, ).")
        if len(args.slot_top_center_x) != 3:
            raise ValueError("slot_top_center_x must be a 3-tuple, (start, stop, num).")
        if len(args.slot_top_center_y) != 3:
            raise ValueError("slot_top_center_y must be a 3-tuple, (start, stop, num).")
        I_fc_pk = args.I_fc[0]
        xstart, xstop, xnum = args.slot_top_center_x
        ystart, ystop, ynum = args.slot_top_center_y
        if xnum != ynum:
            raise ValueError(
                f"Inconsistent slot coordinates: {args.slot_top_center_x}, "
                f"{args.slot_top_center_y}."
            )
        xs = np.linspace(xstart, xstop, int(xnum))
        ys = np.linspace(ystart, ystop, int(ynum))
        slot_coords = np.array([xs, ys]).T
        slot_top_center = slot_coords[args.index]

    npoints = max(1, int(args.points_per_cycle * args.cycles))
    thetas = np.linspace(0, 2 * np.pi * args.cycles, npoints)
    I_fc = I_fc_pk * np.cos(thetas)

    if args.slot_size is None:
        device = make_film(
            args.film_radius,
            args.xi,
            args.lam,
            args.d,
            args.film_points,
            args.film_smooth,
            max_edge_length=args.max_edge_length,
            shape=args.film_shape,
        )
    else:
        device = make_film_with_slot(
            args.film_radius,
            args.xi,
            args.lam,
            args.d,
            args.film_points,
            args.film_smooth,
            args.slot_size,
            slot_top_center,
            args.slot_radius,
            max_edge_length=args.max_edge_length,
            shape=args.film_shape,
        )
    device.layer.gamma = args.gamma
    print(repr(device))

    if args.squid_type == "loop":
        logger.info("Simulating dynamics with 1D loop field coil.")
        if args.pl_center is None:
            pl_center = args.fc_center
        else:
            pl_center = args.pl_center

        simulate_dynamics_loop(
            directory=args.directory,
            device=device,
            fc_radius=args.fc_radius,
            fc_center=args.fc_center,
            pl_radius=args.pl_radius,
            pl_center=pl_center,
            field_units=args.field_units,
            current_units=args.current_units,
            I_fc=I_fc,
            seed_solutions=args.seed_solutions,
            solve_time=args.solve_time,
            dt_init=args.dt_init,
            save_every=args.save_every,
            screening=args.screening,
            metadata=args_as_dict,
        )
    else:
        logger.info(f"Simulating dynamics with {args.squid_type} SQUID.")
        squid = make_squid(
            args.squid_type, args.squid_points, args.squid_smooth, args.squid_angle
        )

        simulate_dynamics_squid(
            directory=args.directory,
            device=device,
            squid=squid,
            squid_position=args.squid_position,
            squid_iterations=args.squid_iterations,
            field_units=args.field_units,
            current_units=args.current_units,
            I_fc=I_fc,
            seed_solutions=args.seed_solutions,
            solve_time=args.solve_time,
            dt_init=args.dt_init,
            save_every=args.save_every,
            screening=args.screening,
            metadata=args_as_dict,
        )


if __name__ == "__main__":
    main()
