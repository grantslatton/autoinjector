"""Generate a dovetail sled and track as STEP solids and STL meshes.

The sled is sized for the BD 1 mL Luer-Lok syringe (309628).

Requires CadQuery 2.8. Run `python autoinjector.py` to write the parts to output/.
All coordinates are in millimeters.
"""

from argparse import ArgumentParser
from collections import defaultdict
from math import atan, cos, hypot, pi, radians, sin, sqrt, tan
from pathlib import Path
from tempfile import TemporaryDirectory

import cadquery as cq
from cadquery.func import loft, vertex
from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCP.IFSelect import IFSelect_RetDone

# X=0 is the plate's outer face, excluding the spikes. Y is centered on the
# original 15 mm block. Z uses the original block's barrel axis.
PLATE_FRONT_X = 0.0
CENTER_Y = 7.5
BARREL_Z = 7.5

PLATE_THICKNESS = 4.0
PLATE_BACK_X = PLATE_FRONT_X + PLATE_THICKNESS
TRACK_LENGTH_FROM_PLATE = 260.0
TRACK_OPEN_END_X = PLATE_FRONT_X + TRACK_LENGTH_FROM_PLATE

# In the seated position, the plate thickness plus the sled's nose-to-slot
# length is the distance from the plate's outside face to the coin slot.
SLED_NOSE_X = PLATE_BACK_X
PLATE_TO_COIN_SLOT = 97.0
COIN_SLOT_FRONT_X = PLATE_FRONT_X + PLATE_TO_COIN_SLOT
COIN_SLOT_WIDTH = 2.5
COIN_SLOT_BACK_X = COIN_SLOT_FRONT_X + COIN_SLOT_WIDTH
SLED_TAIL_LENGTH = 12.5
SLED_BACK_X = COIN_SLOT_BACK_X + SLED_TAIL_LENGTH

PLATE_TO_RELIEF = 88.0
SLED_RELIEF_X = PLATE_FRONT_X + PLATE_TO_RELIEF

WALL_HALF_WIDTH = 8.5
LEFT_WALL, RIGHT_WALL = CENTER_Y - WALL_HALF_WIDTH, CENTER_Y + WALL_HALF_WIDTH
RAIL_REACH = 3.0
LEFT_RAIL, RIGHT_RAIL = LEFT_WALL - RAIL_REACH, RIGHT_WALL + RAIL_REACH
SLED_BOTTOM_Z = -3.5
SLED_RAIL_HEIGHT = 4.5
SLED_RAIL_TOP_Z = SLED_BOTTOM_Z + SLED_RAIL_HEIGHT
SLED_TOP_Z = 16.0

# Size the female dovetail from the male rail. The plate spans the full housing
# width and extends equally above and below the barrel axis.
RAIL_CLEARANCE = 0.30
TRACK_LIP_HEIGHT = 1.5
TRACK_RAIL_TOP_Z = SLED_RAIL_TOP_Z + TRACK_LIP_HEIGHT
TRACK_SIDE_THICKNESS = 5.5
TRACK_HALF_WIDTH = WALL_HALF_WIDTH + RAIL_REACH + RAIL_CLEARANCE + TRACK_SIDE_THICKNESS
TRACK_LEFT_Y, TRACK_RIGHT_Y = CENTER_Y - TRACK_HALF_WIDTH, CENTER_Y + TRACK_HALF_WIDTH
TRACK_BASE_THICKNESS = 4.4
TRACK_BASE_Z = SLED_BOTTOM_Z - TRACK_BASE_THICKNESS
TRACK_TOP_Z = 2 * BARREL_Z - TRACK_BASE_Z

REAR_TAB_LENGTH = 3.0
REAR_TAB_CLEARANCE = 1.0


# Shared cross-section

TOP_ROUND = 2.0


def yz(x: float) -> cq.Workplane:
    """Place a cross-section at station X so it extrudes toward the rear."""
    return cq.Workplane("YZ", origin=(x, 0, 0))


def outside_section(x: float, end_x: float, crown: float = 0.0) -> cq.Solid:
    """Extrude the dovetail and outer walls, with rounded top corners if needed."""
    outline = (
        yz(x)
        .moveTo(LEFT_RAIL, SLED_BOTTOM_Z)
        .lineTo(RIGHT_RAIL, SLED_BOTTOM_Z)
        .lineTo(RIGHT_WALL, SLED_RAIL_TOP_Z)
        .lineTo(RIGHT_WALL, SLED_TOP_Z - crown)
    )
    if crown:
        inset = crown * (1 - 1 / sqrt(2))
        outline = (
            outline.threePointArc(
                (RIGHT_WALL - inset, SLED_TOP_Z - inset),
                (RIGHT_WALL - crown, SLED_TOP_Z),
            )
            .lineTo(LEFT_WALL + crown, SLED_TOP_Z)
            .threePointArc(
                (LEFT_WALL + inset, SLED_TOP_Z - inset), (LEFT_WALL, SLED_TOP_Z - crown)
            )
        )
    else:
        outline = outline.lineTo(LEFT_WALL, SLED_TOP_Z)
    return outline.lineTo(LEFT_WALL, SLED_RAIL_TOP_Z).close().extrude(end_x - x).val()


def barrel(
    x: float, end_x: float, radius: float, *, top: float = SLED_TOP_Z + 1
) -> cq.Solid:
    """Cut the open channel with one semicircle and two upright walls."""
    left, right = CENTER_Y - radius, CENTER_Y + radius
    return (
        yz(x)
        .moveTo(left, top)
        .lineTo(left, BARREL_Z)
        .threePointArc((CENTER_Y, BARREL_Z - radius), (right, BARREL_Z))
        .lineTo(right, top)
        .close()
        .extrude(end_x - x)
        .val()
    )


# Sled nose and spring pockets

CHANNEL_RADIUS = 5.0
NOSE_ROUND = 1.0
SPRING_PAIR_COUNT = 3
POCKET_LENGTH = 10.8
SPRING_END_GAP = 0.5
OUTER_WALL_THICKNESS = 1.0
POCKET_DEPTH_BELOW_BARREL_FLOOR = 0.3
POCKET_ROOF_THICKNESS = 3.0
POCKET_OVERCUT = 3.0
SPRING_ROOT_OVERLAP = 0.2
SPRING_FLOOR_OVERLAP = 0.05


def build_sled_body() -> cq.Shape:
    """Make the rounded front and evenly space the spring pairs before the relief."""
    (stock,) = (
        outside_section(SLED_NOSE_X, SLED_RELIEF_X, crown=TOP_ROUND)
        .cut(barrel(SLED_NOSE_X, SLED_RELIEF_X, CHANNEL_RADIUS))
        .Solids()
    )
    stock = stock.fillet(NOSE_ROUND, stock.faces("<X").Edges())

    # The first and last spring centers are one pitch from the layout ends.
    pitch = (SLED_RELIEF_X - PLATE_FRONT_X) / (SPRING_PAIR_COUNT + 1)
    first_pocket = PLATE_FRONT_X + pitch - POCKET_LENGTH / 2
    wall_y = LEFT_WALL + OUTER_WALL_THICKNESS
    channel_y = CENTER_Y - CHANNEL_RADIUS
    floor_z = BARREL_Z - CHANNEL_RADIUS
    throat = channel_y - wall_y
    pocket_bottom = -POCKET_DEPTH_BELOW_BARREL_FLOOR
    pocket_top = SLED_TOP_Z - POCKET_ROOF_THICKNESS - floor_z

    # The pocket roof follows a 45-degree slope between the wall and channel.
    pocket = (
        cq.Workplane("YZ", origin=(first_pocket, wall_y, floor_z))
        .polyline(
            [
                (0, pocket_bottom),
                (throat + POCKET_OVERCUT, pocket_bottom),
                (throat + POCKET_OVERCUT, pocket_top),
                (throat, pocket_top),
                (0, pocket_top - throat),
            ]
        )
        .close()
        .extrude(POCKET_LENGTH)
        .val()
    )

    # The spring back is measured from the pocket wall. Measure the contact
    # face from the channel wall so a radius change carries the bite with it.
    root_bottom = pocket_bottom - SPRING_FLOOR_OVERLAP
    # Measured profile stations define the stem and contact curve in the YZ sketch.
    stem_y = 1.24
    crest_z = 7.54
    spring = (
        cq.Workplane("YZ", origin=(first_pocket + SPRING_END_GAP, wall_y, floor_z))
        .moveTo(-SPRING_ROOT_OVERLAP, root_bottom)
        .lineTo(-SPRING_ROOT_OVERLAP, 0)
        .lineTo(0, 0)
        .lineTo(0.5689, 0.1422)
        .spline(
            [(0.64, 0.20), (0.90, 0.62), (1.09, 1.12), (1.20, 1.78), (stem_y, 2.55)],
            includeCurrent=True,
        )
        .lineTo(stem_y, 7.37)
        .lineTo(1.43, crest_z)
        .lineTo(throat - 0.38, crest_z)
        .lineTo(throat - 0.1302, 7.42)
        .lineTo(throat + 0.63, 6.68)
        .spline(
            [
                (throat + dy, z)
                for dy, z in (
                    (0.78, 6.38),
                    (0.83, 6.21),
                    (0.84, 6.08),
                    (0.82, 5.95),
                    (0.72, 5.70),
                    (0.4726, 5.41),
                    (0.2963, 5.15),
                )
            ],
            includeCurrent=True,
        )
        .lineTo(throat, 4.35)
        .lineTo(throat, root_bottom)
        .close()
        .extrude(POCKET_LENGTH - 2 * SPRING_END_GAP)
        .val()
    )

    # Cut the mirrored pockets before joining the leaves so the end gaps remain.
    def paired_pattern(shape: cq.Shape) -> list[cq.Shape]:
        """Mirror the cutter or leaf, then place each pair at the spring pitch."""
        pair = (shape, shape.mirror("XZ", (0, CENTER_Y, 0)))
        return [
            side.translate((index * pitch, 0, 0))
            for index in range(SPRING_PAIR_COUNT)
            for side in pair
        ]

    return stock.cut(*paired_pattern(pocket)).fuse(*paired_pattern(spring)).clean()


# Sled coin slot and rear tabs

# The relief clears the flange. The rear barrel is sized separately.
WIDENED_CHANNEL_RADIUS = 5.75
REAR_CHANNEL_RADIUS = 5.5
COIN_LIP_RADIUS = 2.0

# The coin seat keeps the original flange depth when the barrel relief changes.
COIN_SEAT_DEPTH_BELOW_BARREL = 6.5
COIN_SLOT_FLOOR_Z = BARREL_Z - COIN_SEAT_DEPTH_BELOW_BARREL
SIDE_SEAT_RISE = 1.5
SIDE_SEAT_FLAT_WIDTH = 0.5
SEAT_BACK_WALL_HEIGHT = 1.0
# The 22.5° rear face opens toward the top while the seat stays vertical.
SLOT_SLOPE = tan(pi / 8)
SIDE_SEAT_Z = COIN_SLOT_FLOOR_Z + SIDE_SEAT_RISE
SLOT_SLOPE_START_Z = SIDE_SEAT_Z + SEAT_BACK_WALL_HEIGHT
# The seat rises at 45° after the flat land.
SIDE_SEAT_WIDTH = SIDE_SEAT_FLAT_WIDTH + SIDE_SEAT_RISE

TAB_START_X = SLED_BACK_X - REAR_TAB_LENGTH
TAB_UNDERSIDE_SLOPE = 1.0
# The tab climbs across the side clearance before it reaches the track deck.
TAB_ROOT_Z = (
    TRACK_RAIL_TOP_Z + REAR_TAB_CLEARANCE - RAIL_CLEARANCE * TAB_UNDERSIDE_SLOPE
)
TAB_TIP_RADIUS = 3.0
FEATURE_EDGE_RADIUS = 0.5
REAR_PERIMETER_RADIUS = 1.0


def coin_slot_cutter() -> cq.Shape:
    """Make a slot that gives the flange room to tilt above its seat."""
    lip_start_x = COIN_SLOT_FRONT_X - COIN_LIP_RADIUS
    lip_center_z = SLED_TOP_Z - COIN_LIP_RADIUS
    tool_top = SLED_TOP_Z + 1
    span = RIGHT_RAIL - LEFT_RAIL + 2
    section = cq.Workplane(
        cq.Plane(
            (0, RIGHT_RAIL + 1, 0),
            xDir=(1, 0, 0),
            normal=(0, -1, 0),
        )
    )
    opening = (
        section.moveTo(lip_start_x, tool_top)
        .lineTo(lip_start_x, SLED_TOP_Z)
        .threePointArc(
            (
                lip_start_x + COIN_LIP_RADIUS / sqrt(2),
                lip_center_z + COIN_LIP_RADIUS / sqrt(2),
            ),
            (COIN_SLOT_FRONT_X, lip_center_z),
        )
        .lineTo(COIN_SLOT_FRONT_X, COIN_SLOT_FLOOR_Z)
        .lineTo(COIN_SLOT_BACK_X, COIN_SLOT_FLOOR_Z)
        .lineTo(COIN_SLOT_BACK_X, SLOT_SLOPE_START_Z)
        .lineTo(
            COIN_SLOT_BACK_X + (tool_top - SLOT_SLOPE_START_Z) * SLOT_SLOPE, tool_top
        )
        .close()
        .extrude(span)
        .val()
    )

    # Limit the full-width cut so it leaves a flat ledge and ramp at each wall.
    seat_clearance = (
        yz(lip_start_x - 1)
        .polyline(
            [
                (LEFT_RAIL, SIDE_SEAT_Z),
                (LEFT_WALL + SIDE_SEAT_FLAT_WIDTH, SIDE_SEAT_Z),
                (LEFT_WALL + SIDE_SEAT_WIDTH, COIN_SLOT_FLOOR_Z),
                (RIGHT_WALL - SIDE_SEAT_WIDTH, COIN_SLOT_FLOOR_Z),
                (RIGHT_WALL - SIDE_SEAT_FLAT_WIDTH, SIDE_SEAT_Z),
                (RIGHT_RAIL, SIDE_SEAT_Z),
                (RIGHT_RAIL, tool_top),
                (LEFT_RAIL, tool_top),
            ]
        )
        .close()
        .extrude(SLED_BACK_X - lip_start_x + 1)
        .val()
    )
    return opening.intersect(seat_clearance)


def butterfly_tab() -> cq.Shape:
    """Make one rear tab with a rounded tip that meets both sides tangentially."""
    center_z = SLED_TOP_Z - TAB_TIP_RADIUS
    slope = TAB_UNDERSIDE_SLOPE
    # Keep one radius of normal distance between the tip center and the underside.
    center_y = (
        RIGHT_WALL + (center_z - TAB_ROOT_Z - TAB_TIP_RADIUS * hypot(1, slope)) / slope
    )
    ramp_angle = atan(slope)

    def tip(angle: float) -> tuple[float, float]:
        """Locate a point on the arc that joins the top to the ramp."""
        return (
            center_y + TAB_TIP_RADIUS * cos(angle),
            center_z + TAB_TIP_RADIUS * sin(angle),
        )

    return (
        yz(TAB_START_X)
        .moveTo(RIGHT_WALL, TAB_ROOT_Z)
        .lineTo(*tip(ramp_angle - pi / 2))
        .threePointArc(tip(ramp_angle / 2), tip(pi / 2))
        .lineTo(RIGHT_WALL, SLED_TOP_Z)
        .close()
        .extrude(REAR_TAB_LENGTH)
        .val()
    )


def floor_relief_cutter() -> cq.Solid:
    """Cut one half of the small pocket at the back of the coin seat."""
    depth, length = 0.5, 0.05
    reach = RIGHT_WALL - SIDE_SEAT_WIDTH - CENTER_Y
    tool_top = COIN_SLOT_FLOOR_Z + depth
    wire = cq.Wire.makePolygon(
        [
            cq.Vector(COIN_SLOT_BACK_X - length, CENTER_Y, COIN_SLOT_FLOOR_Z),
            cq.Vector(COIN_SLOT_BACK_X, CENTER_Y, COIN_SLOT_FLOOR_Z - depth),
            cq.Vector(COIN_SLOT_BACK_X, CENTER_Y, tool_top),
            cq.Vector(COIN_SLOT_BACK_X - length, CENTER_Y, tool_top),
        ],
        close=True,
    )
    # The oblique extrusion rises to the floor height at the end of the pocket.
    return cq.Solid.extrudeLinear(wire, [], cq.Vector(0, reach, depth))


def fillet_rear_edges(shape: cq.Shape) -> cq.Solid:
    """Round the rear perimeter, coin face and tab roots in one pass."""
    tol = 1e-6
    slope_top = COIN_SLOT_BACK_X + (SLED_TOP_Z - SLOT_SLOPE_START_Z) * SLOT_SLOPE
    wall_y = (LEFT_WALL, RIGHT_WALL)
    channel_y = (CENTER_Y - REAR_CHANNEL_RADIUS, CENTER_Y + REAR_CHANNEL_RADIUS)
    fillet = BRepFilletAPI_MakeFillet(shape.wrapped)
    feature_edges: list[cq.Edge] = []
    perimeter_edges: list[cq.Edge] = []

    def on_plane(extent: tuple[float, float], station: float) -> bool:
        """Both ends of the edge's extent must lie on the same plane."""
        low, high = extent
        return abs(low - station) < tol and abs(high - station) < tol

    for edge in shape.Edges():
        bounds = edge.BoundingBox()
        x_extent = (bounds.xmin, bounds.xmax)
        y_extent = (bounds.ymin, bounds.ymax)

        at_top = bounds.zmin > SLED_TOP_Z - tol
        outboard = bounds.ymin < LEFT_WALL - tol or bounds.ymax > RIGHT_WALL + tol
        if on_plane(x_extent, SLED_BACK_X):
            # Keep the tab rim round small enough to follow its narrow outline.
            tab_rim = bounds.zmin >= TAB_ROOT_Z - tol and (outboard or at_top)
            (feature_edges if tab_rim else perimeter_edges).append(edge)
        else:
            tab_front = on_plane(x_extent, TAB_START_X) and outboard
            on_wall = any(on_plane(y_extent, y) for y in wall_y)
            on_channel = any(on_plane(y_extent, y) for y in channel_y)
            cap_edge = (
                at_top
                and bounds.xmin >= slope_top - tol
                and (on_wall or on_channel or on_plane(x_extent, slope_top))
            )
            seat_edge = (
                on_wall
                and bounds.xmin >= COIN_SLOT_BACK_X - tol
                and bounds.xmax <= slope_top + tol
                and bounds.zmin >= SIDE_SEAT_Z - tol
            )
            if any((tab_front, cap_edge, seat_edge)):
                feature_edges.append(edge)

    counts = (len(feature_edges), len(perimeter_edges))
    if counts != (22, 8):
        raise RuntimeError(
            f"Rear fillet selected {counts[0]} feature edges and {counts[1]} perimeter edges. "
            "Expected 22 feature edges and 8 perimeter edges. Check the slot and tab dimensions."
        )
    for radius, edges in (
        (FEATURE_EDGE_RADIUS, feature_edges),
        (REAR_PERIMETER_RADIUS, perimeter_edges),
    ):
        for edge in edges:
            fillet.Add(radius, edge.wrapped)
    fillet.Build()
    if not fillet.IsDone():
        raise RuntimeError(
            f"OpenCascade could not round the rear sled edges. "
            f"{fillet.NbFaultyContours()} of {fillet.NbContours()} fillet contours failed."
        )
    rounded = cq.Shape.cast(fillet.Shape())
    solids = rounded.Solids()
    if not rounded.isValid() or len(solids) != 1:
        raise RuntimeError("The rear fillets must leave one valid solid.")
    return solids[0]


def build_sled_rear() -> cq.Solid:
    """Cut the rear features from one blank and blend the exposed edges."""
    symmetry = (0, CENTER_Y, 0)
    tab = butterfly_tab()
    # The front keeps the sled's long top round. The tail needs its own edge radii.
    blank = outside_section(SLED_RELIEF_X, COIN_SLOT_FRONT_X, crown=TOP_ROUND).fuse(
        outside_section(COIN_SLOT_FRONT_X, SLED_BACK_X),
        tab,
        tab.mirror("XZ", symmetry),
    )
    relief_cutter = floor_relief_cutter()
    machined = blank.cut(
        barrel(SLED_RELIEF_X, COIN_SLOT_FRONT_X, WIDENED_CHANNEL_RADIUS),
        coin_slot_cutter(),
        barrel(COIN_SLOT_BACK_X, SLED_BACK_X, REAR_CHANNEL_RADIUS),
        relief_cutter,
        relief_cutter.mirror("XZ", symmetry),
    ).clean()
    return fillet_rear_edges(machined)


# Sled wings

WING_TO_COIN_SLOT_GAP = 3.5
WING_WALL_CLEARANCE_ABOVE_DECK = 0.5


def build_wings() -> tuple[cq.Solid, cq.Solid]:
    """Build one wing in wall-local coordinates and mirror it to the other side."""
    # Fillet and trim in a fixed local frame. Moving the coin slot then
    # changes only the placement of the finished wing.
    pad_start, front, end = 0.0, 1.6917, 40.5
    full_width_x = end - 18.0
    pad_width, reach = 1.0, 7.0
    tool_margin = 1.0
    rise, tip_height = 8.0, 1.7
    cross_slope = rise / reach
    # Give the underside cutter room to pass beneath the inboard pad.
    root = (pad_width + tool_margin) * cross_slope + tool_margin
    top = root + rise
    slot_radius, slot_depth, slot_y = 1.5, 3.27, 2.8
    mouth_radius, rim_radius = 0.5, 0.5
    slot_x = end - slot_depth + slot_radius
    near, far = slot_y - slot_radius, slot_y + slot_radius
    mouth_x = end - mouth_radius
    mouth_near, mouth_far = near - mouth_radius, far + mouth_radius
    mid = mouth_radius / sqrt(2)

    # Start with the complete top profile. The pad carries the fillets past the
    # wall during construction so the visible rounds end cleanly at the joint.
    wing = (
        cq.Workplane("XY")
        .moveTo(pad_start, -pad_width)
        .lineTo(end, -pad_width)
        .lineTo(end, mouth_near)
        .threePointArc((mouth_x + mid, mouth_near + mid), (mouth_x, near))
        .lineTo(slot_x, near)
        .threePointArc((slot_x - slot_radius, slot_y), (slot_x, far))
        .lineTo(mouth_x, far)
        .threePointArc((mouth_x + mid, mouth_far - mid), (end, mouth_far))
        .lineTo(end, reach)
        .lineTo(full_width_x, reach)
        .bezier(
            [(full_width_x - 9.5, reach), (front + 5.3083, 0), (front, 0)],
            includeCurrent=True,
        )
        .lineTo(pad_start, 0)
        .close()
        .extrude(top)
        .val()
    )

    # Each underside is a plane. Draw its cutter in the view where its slope is
    # defined. The tip cutter starts below the stock and leaves a short end land.
    y_min, y_max = -pad_width - tool_margin, reach + tool_margin
    x_start, x_end = pad_start - tool_margin, end + tool_margin
    cross_cut = (
        yz(x_start)
        .polyline(
            [
                (y_min, -tool_margin),
                (y_max, -tool_margin),
                (y_max, root + cross_slope * y_max),
                (y_min, root + cross_slope * y_min),
            ]
        )
        .close()
        .extrude(x_end - x_start)
        .val()
    )
    # The cut starts one rise back from the tip and leaves a flat tip face.
    tip_slope = (rise - tip_height) / rise
    tip_z = top - tip_height
    tip_cut = (
        cq.Workplane("XZ", origin=(0, y_max, 0))
        .polyline(
            [
                (end - (tip_z + tool_margin) / tip_slope, -tool_margin),
                (x_end, -tool_margin),
                (x_end, tip_z + tip_slope * tool_margin),
            ]
        )
        .close()
        .extrude(y_max - y_min)
        .val()
    )
    (wing,) = wing.cut(cross_cut, tip_cut).clean().Solids()

    # Find shared edges once. Smooth joins in the outline and edges confined
    # to the disposable pad need no fillet.
    faces_at_edge = defaultdict(list)
    for face in wing.Faces():
        for edge in face.Edges():
            faces_at_edge[edge].append(face)

    rim = []
    for edge in wing.Edges():
        faces = faces_at_edge[edge]
        if len(faces) == 2 and edge.BoundingBox().ymax > 1e-5:
            midpoint = edge.positionAt(0.5)
            a, b = (face.normalAt(midpoint) for face in faces)
            if a.getAngle(b) > radians(2):
                rim.append(edge)

    rounded = wing.fillet(rim_radius, rim)
    trim = cq.Solid.makeBox(x_end - x_start, y_max, top + tool_margin, (x_start, 0, 0))
    exposed = rounded.intersect(trim).Solids()
    if len(exposed) != 1 or not exposed[0].isValid():
        raise ValueError("Rounding and trimming the wing did not leave one valid solid")
    attach_x = COIN_SLOT_FRONT_X - WING_TO_COIN_SLOT_GAP - end
    attach_z = TRACK_RAIL_TOP_Z + WING_WALL_CLEARANCE_ABOVE_DECK - root
    positive = exposed[0].translate((attach_x, RIGHT_WALL, attach_z))
    return positive, positive.mirror("XZ", (0, CENTER_Y, 0))


# Track housing, plate, and fins

RAIL_TOP_ROUND = 1.0
RAIL_BOTTOM_ROUND = 2.0
RAIL_END_ROUND = 1.0

# Offset the groove's two flanks sideways to match the male dovetail.
GROOVE_LEFT_FLOOR, GROOVE_RIGHT_FLOOR = (
    LEFT_RAIL - RAIL_CLEARANCE,
    RIGHT_RAIL + RAIL_CLEARANCE,
)
GROOVE_LEFT_LIP, GROOVE_RIGHT_LIP = (
    LEFT_WALL - RAIL_CLEARANCE,
    RIGHT_WALL + RAIL_CLEARANCE,
)

PLATE_FRONT_ROUND = 2.0
PLATE_REAR_ROUND = 1.0
PLATE_SLOT_DIAMETER = 8.5

SPIKE_PROJECTION = 5.0
SPIKE_HALF_WIDTH = 3.25
SPIKE_HEIGHT = 3.5
SPIKE_PITCH = 8.5

FIN_SETBACK = 15.0
FIN_LENGTH = 40.0
FIN_WIDTH = 4.0
FIN_HEIGHT = 10.0
FIN_EMBED_DEPTH = 2.5
FIN_PLAN_ROUND = 1.0
FIN_CROWN_ROUND = 1.0
FIN_FRONT_TOP_ROUND = 0.5
FIN_SLOT_CENTER_ABOVE_RAIL = 6.5
FIN_SLOT_RADIUS = 1.5
FIN_SLOT_STRAIGHT_LENGTH = 1.0
FIN_SLOT_LIP_ROUND = 0.5

BUTTRESS_LENGTH = 8.0
BUTTRESS_WIDTH = 4.0
BUTTRESS_PLATE_OVERLAP = 0.4
BUTTRESS_RAIL_OVERLAP = 1.5
BUTTRESS_OUTSIDE_ROUND = 1.0

# Measure X from the plate back and Z from the rail top. The first and last
# samples also give the height and length of the reference profile.
BUTTRESS_PROFILE = (
    (0.0, 20.4),
    (0.02, 19.028),
    (0.05, 18.212),
    (0.1, 17.294),
    (0.23, 15.668),
    (0.5, 13.404),
    (1, 10.519),
    (2, 6.618),
    (3, 3.941),
    (4, 2.110),
    (5, 1.008),
    (6, 0.365),
    (7, 0.058),
    (7.5, 0.008),
    (8, 0.0),
)

# Extend tools beyond the faces they act on so the cut runs through the part.
CUT_MARGIN = 0.1
EDGE_LOCATION_TOLERANCE = 1e-6


def rail_housing(open_end_x: float) -> cq.Shape:
    """Extrude the dovetail groove and round the insertion end."""
    # Start where the plate's front round ends to keep its outside edge intact.
    start_x = PLATE_FRONT_X + PLATE_FRONT_ROUND
    top_tangent_z = TRACK_RAIL_TOP_Z - RAIL_TOP_ROUND
    bottom_tangent_z = TRACK_BASE_Z + RAIL_BOTTOM_ROUND
    profile = (
        yz(start_x)
        .moveTo(GROOVE_LEFT_FLOOR, SLED_BOTTOM_Z)
        .lineTo(GROOVE_LEFT_LIP, SLED_RAIL_TOP_Z)
        .lineTo(GROOVE_LEFT_LIP, TRACK_RAIL_TOP_Z)
        .lineTo(TRACK_LEFT_Y + RAIL_TOP_ROUND, TRACK_RAIL_TOP_Z)
        .tangentArcPoint((TRACK_LEFT_Y, top_tangent_z), relative=False)
        .lineTo(TRACK_LEFT_Y, bottom_tangent_z)
        .tangentArcPoint(
            (TRACK_LEFT_Y + RAIL_BOTTOM_ROUND, TRACK_BASE_Z), relative=False
        )
        .lineTo(TRACK_RIGHT_Y - RAIL_BOTTOM_ROUND, TRACK_BASE_Z)
        .tangentArcPoint((TRACK_RIGHT_Y, bottom_tangent_z), relative=False)
        .lineTo(TRACK_RIGHT_Y, top_tangent_z)
        .tangentArcPoint(
            (TRACK_RIGHT_Y - RAIL_TOP_ROUND, TRACK_RAIL_TOP_Z), relative=False
        )
        .lineTo(GROOVE_RIGHT_LIP, TRACK_RAIL_TOP_Z)
        .lineTo(GROOVE_RIGHT_LIP, SLED_RAIL_TOP_Z)
        .lineTo(GROOVE_RIGHT_FLOOR, SLED_BOTTOM_Z)
        .close()
    )
    return (
        profile.extrude(open_end_x - start_x).edges(">X").fillet(RAIL_END_ROUND).val()
    )


def slotted_plate() -> cq.Shape:
    """Build the end plate with a U-slot centered on the sled barrel."""
    block = cq.Solid.makeBox(
        PLATE_THICKNESS,
        TRACK_RIGHT_Y - TRACK_LEFT_Y,
        TRACK_TOP_Z - TRACK_BASE_Z,
        cq.Vector(PLATE_FRONT_X, TRACK_LEFT_Y, TRACK_BASE_Z),
    )
    # Round the front and its lengthwise corners together, then the smaller rear rim.
    front = block.fillet(PLATE_FRONT_ROUND, block.edges("<X or |X"))
    rounded = front.fillet(PLATE_REAR_ROUND, front.edges(">X"))

    slot = barrel(
        PLATE_FRONT_X - CUT_MARGIN,
        PLATE_BACK_X + CUT_MARGIN,
        radius=PLATE_SLOT_DIAMETER / 2,
        top=TRACK_TOP_Z + CUT_MARGIN,
    )
    return rounded.cut(slot)


def half_diamond_spike() -> cq.Shape:
    """Loft a triangular root to a point along negative X."""
    root = (
        cq.Workplane("YZ")
        .polyline([(-SPIKE_HALF_WIDTH, 0), (0, -SPIKE_HEIGHT), (SPIKE_HALF_WIDTH, 0)])
        .close()
        .val()
    )
    return loft(root, vertex(-SPIKE_PROJECTION, 0, 0), cap=True, ruled=True)


def rail_fin(center_y: float) -> cq.Shape:
    """Build a rounded fin with a sideways U-slot at its raised end."""
    outer_y = center_y + FIN_WIDTH / 2
    front_x = PLATE_BACK_X + FIN_SETBACK
    rear_x = front_x + FIN_LENGTH
    top_z = TRACK_RAIL_TOP_Z + FIN_HEIGHT
    root_z = TRACK_RAIL_TOP_Z - FIN_EMBED_DEPTH
    xz = cq.Plane((0, outer_y, 0), (1, 0, 0), (0, -1, 0))

    # Three control points at each height give the quintic roof zero slope
    # and zero curvature where it starts and where it meets the rail.
    heights = (
        top_z,
        top_z,
        top_z,
        TRACK_RAIL_TOP_Z,
        TRACK_RAIL_TOP_Z,
        TRACK_RAIL_TOP_Z,
    )
    poles = [(front_x + FIN_LENGTH * i / 5, z) for i, z in enumerate(heights)]
    blank = (
        cq.Workplane(xz)
        .moveTo(front_x, root_z)
        .lineTo(*poles[0])
        .bezier(poles[1:], includeCurrent=True)
        .lineTo(rear_x, root_z)
        .close()
        .extrude(FIN_WIDTH)
        .val()
    )
    roof_edges = [edge for edge in blank.Edges() if edge.geomType() == "BEZIER"]
    rounded = blank.fillet(FIN_CROWN_ROUND, roof_edges)
    front_top_edges = [
        edge
        for edge in rounded.Edges()
        if abs(edge.Center().x - front_x) < EDGE_LOCATION_TOLERANCE
        and abs(edge.Center().z - top_z) < EDGE_LOCATION_TOLERANCE
    ]
    rounded = rounded.fillet(FIN_FRONT_TOP_ROUND, front_top_edges)

    # Round the front corners in plan to keep their axes vertical at the roof.
    nose_tangent_x = front_x + FIN_PLAN_ROUND
    inner_y = outer_y - FIN_WIDTH
    footprint = (
        cq.Workplane("XY", origin=(0, 0, root_z - CUT_MARGIN))
        .moveTo(nose_tangent_x, outer_y)
        .lineTo(rear_x + CUT_MARGIN, outer_y)
        .lineTo(rear_x + CUT_MARGIN, inner_y)
        .lineTo(nose_tangent_x, inner_y)
        .tangentArcPoint((front_x, inner_y + FIN_PLAN_ROUND), relative=False)
        .lineTo(front_x, outer_y - FIN_PLAN_ROUND)
        .tangentArcPoint((nose_tangent_x, outer_y), relative=False)
        .close()
        .extrude(top_z - root_z + 2 * CUT_MARGIN)
        .val()
    )
    rounded = rounded.intersect(footprint)
    lip_radius = FIN_SLOT_LIP_ROUND
    lip_x = front_x + lip_radius
    slot_x = front_x + FIN_SLOT_STRAIGHT_LENGTH
    slot_z = TRACK_RAIL_TOP_Z + FIN_SLOT_CENTER_ABOVE_RAIL
    slot_low_z, slot_high_z = slot_z - FIN_SLOT_RADIUS, slot_z + FIN_SLOT_RADIUS

    # Put the entrance rounds in the slot sketch to preserve the semicircle and
    # straight walls behind them.
    slot = (
        cq.Workplane(xz)
        .workplane(offset=-CUT_MARGIN)
        .moveTo(front_x - CUT_MARGIN, slot_low_z - lip_radius)
        .lineTo(front_x, slot_low_z - lip_radius)
        .radiusArc((lip_x, slot_low_z), lip_radius)
        .lineTo(slot_x, slot_low_z)
        .tangentArcPoint((slot_x, slot_high_z), relative=False)
        .lineTo(lip_x, slot_high_z)
        .radiusArc((front_x, slot_high_z + lip_radius), lip_radius)
        .lineTo(front_x - CUT_MARGIN, slot_high_z + lip_radius)
        .close()
        .extrude(FIN_WIDTH + 2 * CUT_MARGIN)
        .val()
    )
    return rounded.cut(slot)


def plate_buttress() -> cq.Shape:
    """Build the right-side brace tangent to the plate and rail top."""
    plane = cq.Plane((0, TRACK_RIGHT_Y, 0), (1, 0, 0), (0, -1, 0))

    # Scale the measured curve to meet both parts when their heights change.
    length_scale = BUTTRESS_LENGTH / BUTTRESS_PROFILE[-1][0]
    height_scale = (TRACK_TOP_Z - TRACK_RAIL_TOP_Z) / BUTTRESS_PROFILE[0][1]
    rise = [
        (PLATE_BACK_X + x * length_scale, TRACK_RAIL_TOP_Z + z * height_scale)
        for x, z in BUTTRESS_PROFILE[1:]
    ]
    heel = PLATE_BACK_X - BUTTRESS_PLATE_OVERLAP
    toe = PLATE_BACK_X + BUTTRESS_LENGTH
    root = TRACK_RAIL_TOP_Z - BUTTRESS_RAIL_OVERLAP
    blank = (
        cq.Workplane(plane)
        .moveTo(heel, root)
        .lineTo(heel, TRACK_TOP_Z)
        .lineTo(PLATE_BACK_X, TRACK_TOP_Z)
        .spline(rise, includeCurrent=True, tangents=((0, -1), (1, 0)))
        .lineTo(toe, root)
        .close()
        .extrude(BUTTRESS_WIDTH)
        .val()
    )
    outer_roof_edges = [
        edge
        for edge in blank.Edges()
        if edge.geomType() == "BSPLINE"
        and edge.BoundingBox().ymin > TRACK_RIGHT_Y - EDGE_LOCATION_TOLERANCE
    ]
    rounded = blank.fillet(BUTTRESS_OUTSIDE_ROUND, outer_roof_edges)

    # Trim the brace below the rear rim so it cannot fill the plate's round.
    clip_bottom = root - CUT_MARGIN
    clip_top = TRACK_TOP_Z - PLATE_REAR_ROUND
    clip = cq.Solid.makeBox(
        toe - heel + 2 * CUT_MARGIN,
        BUTTRESS_WIDTH + 2 * CUT_MARGIN,
        clip_top - clip_bottom,
        cq.Vector(
            heel - CUT_MARGIN, TRACK_RIGHT_Y - BUTTRESS_WIDTH - CUT_MARGIN, clip_bottom
        ),
    )
    return rounded.intersect(clip)


def build_track(open_end_x: float = TRACK_OPEN_END_X) -> cq.Shape:
    """Fuse the housing, end plate, and repeated features into one track."""
    right_buttress = plate_buttress()

    # Center each fin on the rail before rounding its outer corner.
    fin_center = (GROOVE_RIGHT_LIP + TRACK_RIGHT_Y) / 2
    right_fin = rail_fin(fin_center)
    fin_pitch = 2 * (fin_center - CENTER_Y)

    spike = half_diamond_spike()
    stations = (-SPIKE_PITCH, 0, SPIKE_PITCH)
    # Keep the middle and upper center stations clear for the open slot.
    spikes = [
        spike.translate((PLATE_FRONT_X, CENTER_Y + y, BARREL_Z + z))
        for y in stations
        for z in stations
        if y != 0 or z < 0
    ]
    return (
        rail_housing(open_end_x)
        .fuse(
            slotted_plate(),
            right_buttress,
            right_buttress.mirror("XZ", (0, CENTER_Y, 0)),
            right_fin,
            right_fin.translate((0, -fin_pitch, 0)),
            *spikes,
        )
        .clean()
    )


# Assembly and export

# Absolute deflection in mm and radians keeps long rails and small rounds
# subject to the same mesh tolerances.
STL_LINEAR_DEFLECTION = 0.008
STL_ANGULAR_DEFLECTION = 0.15


def solid(shape: cq.Shape | cq.Workplane, name: str) -> cq.Solid:
    """Require one valid solid so an export cannot contain detached features."""
    if isinstance(shape, cq.Workplane):
        shape = shape.val()
    solids = shape.Solids()
    if len(solids) != 1:
        raise ValueError(
            f"{name} must be one solid, but the build produced {len(solids)}"
        )
    if not solids[0].isValid():
        raise ValueError(f"{name} has invalid CAD geometry")
    return solids[0]


def build_sled() -> cq.Solid:
    """Join the front, rear and both wings in the assembly frame."""
    body = solid(build_sled_body(), "sled body")
    rear = solid(build_sled_rear(), "sled rear")
    wings = build_wings()
    return solid(body.fuse(rear, *wings).clean(), "sled")


def main() -> None:
    """Write both parts only after the full set of exports succeeds."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "output",
        help="destination for the STEP and STL files (default: output/ beside this script)",
    )
    args = parser.parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    parts = {"sled": build_sled(), "track": solid(build_track(), "track")}

    # Stage the files on the destination filesystem. A failed export leaves
    # the previous files available, and a file becomes visible only once complete.
    with TemporaryDirectory(prefix=".build-", dir=output_dir) as staging:
        stage = Path(staging)
        for name, part in parts.items():
            step = stage / f"{name}.step"
            if part.exportStep(str(step)) != IFSelect_RetDone:
                raise OSError(f"Could not export {name} as STEP")
            stl = stage / f"{name}.stl"
            if not part.exportStl(
                str(stl),
                tolerance=STL_LINEAR_DEFLECTION,
                angularTolerance=STL_ANGULAR_DEFLECTION,
                relative=False,
            ):
                raise OSError(f"Could not export {name} as STL")
        for path in sorted(stage.iterdir()):
            path.replace(output_dir / path.name)

    for name, part in parts.items():
        print(f"{name}: one valid solid, {part.Volume():.2f} mm³", flush=True)


if __name__ == "__main__":
    main()
