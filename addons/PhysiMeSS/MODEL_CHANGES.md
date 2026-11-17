# PhysiMeSS: changes of the model since version 1.0.1

This file lists every change of behaviour of PhysiMeSS relative to the published version
(1.0.1, [Gigabyte 2024](https://doi.org/10.46471/gigabyte.136)): what changed, why, and what it
does to the results. Changes which only reorganise the code are not listed.

Simulations made with this version are **not** expected to reproduce the numbers of version 1.0.1
when cells push or rotate fibres. Fibre initialisation in 2D, cell-fibre adhesion and fibre
degradation are unchanged unless stated below.

The geometry used below is in `PhysiMeSS_geometry.h` and tested in `tests/unit/physimess`
(`make physimess-tests`).

## Notation

A fibre is a segment of length `L` with centre `c` and unit orientation `u`. A cell acts on a fibre
at the point `p` of the fibre axis nearest to the cell centre. The fibre rotates about a pivot `o`:
its centre if it has no crosslink, its crosslink point if it has exactly one. `r = p - o` is the
lever arm and `d = |c - o|` the distance from the pivot to the centre. Fibres with two or more
crosslinks do not move (unchanged).

## 1. Fibres in 3D

*Contributed by Ambre Jaeger.*

- **What.** With `use_2D` false, fibres are oriented, pushed and rotated in three dimensions.
  Anisotropic fibres take their orientation from two angles: `fibre_angle` (in the xy plane) and
  `fibre_angle_phi` (from the z axis). Crosslinks between fibres are detected in 3D.
- **What.** Fibres can be created at random until they fill a given fraction of the domain
  (`<fibres_from_file>` with `<enable>false</enable>` and `<relative_fibre_volume>`), instead of
  being read from the csv file.
- **What.** Fibres which overlap a cell at initialisation are removed, in 2D and in 3D.
- **Effect.** 2D simulations in which fibres overlapped cells at time zero start with fewer
  fibres. The hinge example was affected (its vertical fibre was removed): its cell now starts at
  x = 18 instead of x = 25.

## 2. Rotation of a fibre by a cell

This is the part of the model which changed most. In version 1.0.1 (2D only), each cell in contact
with a fibre rotated it immediately by the angle

    angle = fibre_sticky * migration_speed * |r| / (L^2 / 2)

always anticlockwise, whatever the side of the fibre and the direction of the cell.

### 2.1 Direction of the rotation

- **What.** The rotation is about the axis of the torque `r x v`, where `v` is the velocity given
  by the cell (2.4): the point of the fibre touched by the cell follows it, in 2D and in 3D.
- **Why.** A rotation which is always anticlockwise is not a mechanical response to the cell.
- **Effect.** Fibres can turn clockwise; in 3D they turn out of the plane when the cell does not
  move in the plane of the fibre.

### 2.2 Angular velocity, one rotation per mechanics step

- **What.** Each cell in contact with a fibre adds an angular velocity `omega` to the fibre. The
  fibre is rotated once per mechanics step, by the angle `|sum of omega| * dt_mechanics` about the
  axis `sum of omega`, with Rodrigues' formula
  `v' = v cos(a) + (k x v) sin(a) + k (k . v)(1 - cos(a))`
  applied to the orientation and to the vector from the pivot to the centre.
- **Why.** The angle did not depend on the time step, so the rotation rate was inversely
  proportional to `dt_mechanics`. Cells rotated the fibre one after the other, so the result
  depended on the order of the cells, and two cells with opposite torques did not cancel exactly.
- **Effect.** Results no longer depend on `dt_mechanics` beyond the usual discretisation error
  (3D rotating example: final fibre angle 108.2, 107.7 and 106.7 degrees for dt = 0.1, 0.05 and
  0.02 min). The rotation is applied at the beginning of the next mechanics step: saved outputs
  lag the rotation by one mechanics step.

### 2.3 Resistance of the fibre to rotation

- **What.** `omega = (r x v) / (L^2/12 + d^2)`.
- **Why.** PhysiCell mechanics is overdamped (velocities are proportional to forces). For a rod
  with a uniform drag along its length, a force which would translate the rod at the velocity `v`,
  applied at the lever arm `r`, rotates it at this angular velocity: `L^2/12 + d^2` is the mean
  square distance of the points of the rod to the pivot. Version 1.0.1 used `L^2/2` and the
  intermediate 3D code `L^2/12`, whatever the pivot.
- **Effect.** A free fibre is 6 times easier to rotate than in version 1.0.1. A fibre hinged at
  its end is 4 times harder to rotate than about its centre.

### 2.4 Force of the cell: contact and traction

- **What.** The velocity `v` given by a cell to a fibre is the sum of
  - the *contact*: the repulsion of PhysiCell, `sqrt(repulsion_cell * repulsion_fibre) (1 - D/R)^2`,
    directed from the cell to the fibre (`D`: distance from the fibre axis to the cell centre,
    `R`: cell radius + fibre radius). It is the velocity already used to push fibres;
  - the *traction*: `fibre_sticky * motility_vector` of the cell: a migrating cell drags the
    fibres it adheres to (2.5).
- **Why.** In version 1.0.1 only the migration of the cell rotated a fibre: a crowd of cells
  growing against a fibre could translate it but never turn it, and the depth of the contact
  played no role. Cells both displace fibres sterically and pull on them through adhesions; the
  two effects now have one parameter each (repulsion strength, `fibre_sticky`).
- **Effect.** Non-motile cells pressing on a fibre rotate it. With `fibre_sticky = 0` only the
  contact acts. In the cell-fibre mechanics examples, where a single migrating cell meets a
  fibre, the traction dominates: in the 3D hinge example the fibre is rotated by about 10 degrees
  by the contact alone (`fibre_sticky = 0`) and by about 60 degrees with `fibre_sticky = 1`.
- **Flags.** `fibre_rotation` (cell custom data) switches the rotation, by contact and traction;
  `fibre_pushing` switches the translation by contact, as before.

### 2.5 Meaning and value of `fibre_sticky`

- **What.** `fibre_sticky` is a custom data of the *cell*: it represents the adhesion of the cell
  to fibres (the level of its adhesion receptors), and multiplies the traction only (2.4). With 0
  the cell does not adhere: it rotates fibres by contact only. With 1 it drags a fibre with the
  force which would translate the fibre at the velocity of the cell. Cell types can have
  different values.
- **Why.** In version 1.0.1 it was already read from the cell, as a factor of the rotation angle
  per step. The intermediate 3D code read it from the fibre, so that all cell types rotated a
  given fibre equally.
- **Effect.** The default stays 1, as in version 1.0.1, but its effect is not the same (2.2, 2.3).
  The outcome of the examples depends weakly on values above 1 (3D hinge: about 60 degrees with
  1, 65 with 10), because the fibre stops rotating when the cell has passed.

### 2.6 Fibres with one crosslink

- **What.** The fibre rotates rigidly about its crosslink point: orientation and centre are
  rotated by the same rotation, so the crosslink stays at the same place on the fibre.
- **Why.** In the intermediate 3D code the centre moved in the wrong direction when the cell
  pushed between the crosslink and the near end of the fibre.
- **Effect.** The hinge examples behave as hinges: the centre of the hinged fibre stays at a
  constant distance from the crosslink.

### 2.7 Orientation vector

- **What.** The orientation of a fibre keeps its sign when the fibre rotates.
- **Effect.** None on the mechanics (a fibre and its opposite are the same segment); the saved
  `orientation` no longer flips.

## 3. Pushing

- **What.** A cell in contact with a fibre pushes it even when it does not rotate it (contact at
  the pivot, or cell moving along the fibre).
- **Why.** In the intermediate 3D code the function returned before the push in these cases.

## 4. Mechanics time step

- **What.** The fibre mechanics (neighbours, crosslinks, rotation) runs every `dt_mechanics`,
  whatever its value.
- **Why.** In version 1.0.1 it only ran when `dt_mechanics` was 0.1 min: with another value cells
  and fibres did not interact at all.
- **Note.** PhysiCell itself requires `dt_mechanics` to be a multiple of `dt_diffusion`.

## 5. Robustness

- A fibre of zero length (the length is drawn from a normal distribution) no longer gives NaN
  for the nearest point: it is treated as a point.
- Fibres overlapping cells at initialisation are collected first and deleted afterwards (they
  were deleted while iterating over the agents).

## Known limitations

- The rotation angle per step is not capped.
- Fibres named `fibre_vertical` and `fibre_horizontal` have a fixed orientation whatever the
  configuration says (as in version 1.0.1).
- 3D configurations with anisotropic fibres must define `fibre_angle_phi`; when it is missing,
  the first custom variable is silently used instead.
- Fibres with one crosslink are rotated but not translated, except pairs of fibres crosslinked
  to each other alone (as in version 1.0.1).
