#ifndef __PhysiMeSS_geometry_h__
#define __PhysiMeSS_geometry_h__

/*
    Pure geometry used by the cell-fibre mechanics of PhysiMeSS.

    Nothing in this header depends on PhysiCell or BioFVM, so the functions
    can be unit tested on their own (see tests/unit/physimess).

    A fibre is a line segment described by its centre, a unit orientation
    vector and its half length: it spans
        centre - half_length * orientation  ...  centre + half_length * orientation
    All vectors have three components.
*/

#include <algorithm>
#include <cmath>
#include <vector>

namespace PhysiMeSS_Geometry
{

inline double dot( const std::vector<double>& a, const std::vector<double>& b )
{
    double out = 0.0;
    for (unsigned int i = 0; i < 3; i++) {
        out += a[i] * b[i];
    }
    return out;
}

inline double norm( const std::vector<double>& a )
{
    return std::sqrt(dot(a, a));
}

inline std::vector<double> cross( const std::vector<double>& a, const std::vector<double>& b )
{
    std::vector<double> out(3, 0.0);
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
    return out;
}

/*
    Unit orientation of a fibre from its two angles, in radians:
    theta, the angle from the x axis of its projection on the xy plane, and
    phi, its angle from the z axis (spherical coordinates):
        (cos(theta) sin(phi), sin(theta) sin(phi), cos(phi))
    With phi = pi/2 the fibre is in the xy plane, at the angle theta used in
    2D simulations.
*/
inline std::vector<double> orientation_from_angles( double theta, double phi )
{
    std::vector<double> orientation(3, 0.0);
    orientation[0] = std::cos(theta) * std::sin(phi);
    orientation[1] = std::sin(theta) * std::sin(phi);
    orientation[2] = std::cos(phi);
    return orientation;
}

/*
    Vector from the point of the fibre nearest to `point`, to `point`
    (so that nearest_point = point - displacement).

    The nearest point is clamped to the segment: beyond either end of the
    fibre, it is the corresponding endpoint.
*/
inline std::vector<double> displacement_from_fibre(
    const std::vector<double>& point,
    const std::vector<double>& centre,
    const std::vector<double>& orientation,
    double half_length )
{
    const double fibre_length = 2 * half_length;

    std::vector<double> displacement(3, 0.0);
    if (half_length <= 0.) {
        // fibre reduced to a point
        for (unsigned int i = 0; i < 3; i++) {
            displacement[i] = point[i] - centre[i];
        }
        return displacement;
    }

    std::vector<double> endpoint(3, 0.0);          // "start" end of the fibre
    std::vector<double> endpoint_to_point(3, 0.0);
    // endpoint_to_point . (end - start): position along the fibre, scaled by the fibre length
    double projection = 0.0;

    for (unsigned int i = 0; i < 3; i++) {
        endpoint[i] = centre[i] - half_length * orientation[i];
        endpoint_to_point[i] = point[i] - endpoint[i];
        projection += endpoint_to_point[i] * 2 * (centre[i] - endpoint[i]);
    }

    if (projection < 0.) {
        // nearest to the start endpoint
        for (unsigned int i = 0; i < 3; i++) {
            displacement[i] = endpoint_to_point[i];
        }
    }
    else if (projection > fibre_length * fibre_length) {
        // nearest to the other endpoint
        for (unsigned int i = 0; i < 3; i++) {
            displacement[i] = point[i] - (centre[i] + half_length * orientation[i]);
        }
    }
    else {
        // nearest to a point along the fibre
        for (unsigned int i = 0; i < 3; i++) {
            displacement[i] = endpoint_to_point[i]
                - (centre[i] - endpoint[i]) / half_length * projection / fibre_length;
        }
    }
    return displacement;
}

/*
    Velocity given to a fibre by the contact with a cell: the repulsion of
    PhysiCell,
        repulsion * (1 - d / contact_distance)^2
    directed from the cell to the fibre, where d = |displacement| is the
    distance from the fibre axis to the cell centre (`displacement` is the
    result of displacement_from_fibre for the cell centre), and
    contact_distance = cell radius + fibre radius.

    Zero when the cell and the fibre do not touch (d >= contact_distance).
*/
inline std::vector<double> contact_velocity(
    const std::vector<double>& displacement,
    double contact_distance,
    double repulsion )
{
    std::vector<double> velocity(3, 0.0);
    const double distance = std::max(norm(displacement), 0.00001);
    if (distance >= contact_distance) {
        return velocity;
    }
    double magnitude = 1 - distance / contact_distance;
    magnitude *= magnitude; // (1 - d/R)^2
    magnitude *= repulsion;
    magnitude /= distance;
    for (unsigned int i = 0; i < 3; i++) {
        velocity[i] = -magnitude * displacement[i];
    }
    return velocity;
}

/*
    Angular velocity given to a fibre by a cell which acts on it at `impact_point`
    (a point on the fibre axis). `velocity` is the velocity that the force of
    the cell would give to the fibre if it only translated it: the velocity
    of the contact (contact_velocity) and/or of the traction of a migrating
    cell. The fibre rotates about `pivot`: its centre if it is free, its
    crosslink point otherwise.

    The mechanics is overdamped, as for the cells of PhysiCell: velocities are
    proportional to forces. For a rod of length L with a uniform drag along its
    length, a force which would translate the rod at the velocity v, applied
    at the lever arm r = impact_point - pivot, rotates it at

        omega = (r x v) / (L^2/12 + d^2)        with d = |centre - pivot|

    L^2/12 + d^2 is the mean square distance of the points of the rod to the
    pivot (L^2/12 about the centre, L^2/3 about an end): the further the pivot
    is from the centre, the harder it is to rotate the fibre.

    The contributions of several cells add up. Only the part of the velocity
    perpendicular to the fibre matters; the result is perpendicular to the fibre.
*/
inline std::vector<double> angular_velocity_about_pivot(
    const std::vector<double>& pivot,
    const std::vector<double>& centre,
    const std::vector<double>& impact_point,
    const std::vector<double>& velocity,
    double half_length )
{
    std::vector<double> lever(3, 0.0);
    std::vector<double> pivot_to_centre(3, 0.0);
    for (unsigned int i = 0; i < 3; i++) {
        lever[i] = impact_point[i] - pivot[i];
        pivot_to_centre[i] = centre[i] - pivot[i];
    }
    const double fibre_length = 2 * half_length;
    const double rotational_drag = fibre_length * fibre_length / 12 + dot(pivot_to_centre, pivot_to_centre);

    std::vector<double> angular_velocity(3, 0.0);
    if (rotational_drag <= 0.) {
        return angular_velocity; // fibre reduced to a point
    }
    angular_velocity = cross(lever, velocity);
    for (unsigned int i = 0; i < 3; i++) {
        angular_velocity[i] /= rotational_drag;
    }
    return angular_velocity;
}

struct FibreRotation
{
    bool rotated = false;             // false: the fibre is left untouched
    std::vector<double> orientation;  // new unit orientation (only if rotated)
    std::vector<double> centre;       // new fibre centre (only if rotated)
};

/*
    Rigid rotation of a fibre about `pivot` by the rotation vector `rotation`:
    its direction is the axis of the rotation and its norm the angle in radians
    (for an angular velocity omega applied during dt, rotation = omega * dt).

    Both the orientation and the vector from the pivot to the centre are
    rotated with Rodrigues' formula

        v_new = v cos(angle) + (k x v) sin(angle) + k (k . v) (1 - cos(angle))

    so the pivot stays at the same place on the fibre, and the orientation
    keeps its sign. The centre only moves when the pivot is not the centre,
    i.e. for a crosslinked fibre.

    No rotation is returned when the angle is negligible.
*/
inline FibreRotation rotate_fibre_about_pivot(
    const std::vector<double>& pivot,
    const std::vector<double>& centre,
    const std::vector<double>& orientation,
    const std::vector<double>& rotation )
{
    FibreRotation result;

    const double angle = norm(rotation);
    if (angle < 1e-16) {
        return result;
    }
    std::vector<double> axis(3, 0.0);
    std::vector<double> pivot_to_centre(3, 0.0);
    for (unsigned int i = 0; i < 3; i++) {
        axis[i] = rotation[i] / angle;
        pivot_to_centre[i] = centre[i] - pivot[i];
    }

    const double c = std::cos(angle);
    const double s = std::sin(angle);

    const std::vector<double> axis_cross_orientation = cross(axis, orientation);
    const double axis_dot_orientation = dot(axis, orientation);
    const std::vector<double> axis_cross_centre = cross(axis, pivot_to_centre);
    const double axis_dot_centre = dot(axis, pivot_to_centre);

    result.rotated = true;
    result.orientation.assign(3, 0.0);
    result.centre.assign(3, 0.0);
    double new_norm = 1e-32; // same guard as BioFVM::normalize
    for (unsigned int i = 0; i < 3; i++) {
        result.orientation[i] = orientation[i] * c + axis_cross_orientation[i] * s
                                + axis[i] * axis_dot_orientation * (1 - c);
        new_norm += result.orientation[i] * result.orientation[i];
        result.centre[i] = pivot[i] + pivot_to_centre[i] * c + axis_cross_centre[i] * s
                           + axis[i] * axis_dot_centre * (1 - c);
    }
    new_norm = std::sqrt(new_norm);
    for (unsigned int i = 0; i < 3; i++) {
        result.orientation[i] /= new_norm;
    }
    return result;
}

/*
    Effect of a fibre on the velocity of a cell which is within its adhesion
    distance. `previous_velocity` is the velocity of the cell at the previous
    mechanics step and xi = |previous_velocity . orientation| / |previous_velocity|
    measures how much the cell moves along the fibre (1) or across it (0).

    Adhesion: velocity added to the cell, along the fibre and in the direction
    in which the cell already moves,
        vel_adhesion * xi * max(0, 1 - |previous_velocity| / velocity_max)
    The cell is accelerated along the fibre, less and less as its speed
    approaches velocity_max, and not at all above it.
*/
inline std::vector<double> fibre_adhesion_velocity(
    const std::vector<double>& previous_velocity,
    const std::vector<double>& orientation,
    double vel_adhesion,
    double velocity_max )
{
    std::vector<double> adhesion(3, 0.0);
    const double along_fibre = dot(previous_velocity, orientation);
    const double speed = std::max(norm(previous_velocity), 1e-8);
    const double xi = std::fabs(along_fibre) / speed;
    const double magnitude = vel_adhesion * xi * std::max(0.0, 1 - speed / velocity_max);
    const double direction = (along_fibre < 0) ? -1.0 : 1.0;
    for (unsigned int i = 0; i < 3; i++) {
        adhesion[i] = direction * magnitude * orientation[i];
    }
    return adhesion;
}

/*
    Contact: coefficient, between 0 and vel_contact, by which a fibre slows down
    a cell moving across it,
        vel_contact * (1 - xi^2)
    (see damp_motility_across_fibre).
*/
inline double fibre_contact_damping(
    const std::vector<double>& previous_velocity,
    const std::vector<double>& orientation,
    double vel_contact )
{
    const double speed = std::max(norm(previous_velocity), 1e-8);
    const double xi = std::fabs(dot(previous_velocity, orientation)) / speed;
    return vel_contact * (1 - xi * xi);
}

/*
    Motility vector of a cell slowed down by a fibre: the part of the motility
    across the fibre is multiplied by (1 - damping), with damping limited to
    [0, 1]; the part along the fibre is kept.
*/
inline std::vector<double> damp_motility_across_fibre(
    const std::vector<double>& motility,
    const std::vector<double>& orientation,
    double damping )
{
    const double c = std::min(1.0, std::max(0.0, damping));
    const double along_fibre = dot(motility, orientation);
    std::vector<double> damped(3, 0.0);
    for (unsigned int i = 0; i < 3; i++) {
        const double parallel = along_fibre * orientation[i];
        damped[i] = parallel + (1 - c) * (motility[i] - parallel);
    }
    return damped;
}

} // namespace PhysiMeSS_Geometry

#endif
