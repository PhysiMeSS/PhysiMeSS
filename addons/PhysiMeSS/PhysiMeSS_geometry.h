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

struct FibreRotation
{
    bool rotated = false;             // false: the fibre is left untouched
    std::vector<double> orientation;  // new unit orientation (only if rotated)
    std::vector<double> centre;       // new fibre centre (only if rotated)
};

/*
    Rotation of a fibre about `pivot` when `force` is applied at `impact_point`
    (a point on the fibre axis).

    torque  = (impact_point - pivot) x force
    angle   = stickiness * |torque| / (L^2 / 12)        with L = 2 * half_length
    The lever arm (impact_point - pivot) is rotated by `angle` about the torque
    direction with Rodrigues' formula. As the torque is perpendicular to the
    lever arm, the formula reduces to
        v_new = v cos(angle) + (k x v) sin(angle)

    The returned orientation is the direction of the rotated lever arm, i.e.
    it points from the pivot towards the impact point. The fibre centre is
    rotated about the pivot with the rest of the fibre (it only moves when
    the pivot is not the centre, i.e. for a crosslinked fibre).

    No rotation is returned when the force is applied at the pivot or when
    the resulting angle is negligible (force parallel to the lever arm).
*/
inline FibreRotation rotate_fibre_about_pivot(
    const std::vector<double>& pivot,
    const std::vector<double>& centre,
    const std::vector<double>& impact_point,
    const std::vector<double>& force,
    double half_length,
    double stickiness )
{
    FibreRotation result;

    std::vector<double> lever(3, 0.0);
    for (unsigned int i = 0; i < 3; i++) {
        lever[i] = impact_point[i] - pivot[i];
    }
    const double lever_length = norm(lever);
    if (lever_length < 1e-10) {
        return result;
    }

    std::vector<double> old_direction(3, 0.0);
    for (unsigned int i = 0; i < 3; i++) {
        old_direction[i] = lever[i] / lever_length;
    }

    std::vector<double> axis = cross(lever, force);
    const double torque_magnitude = norm(axis);
    const double fibre_length = 2 * half_length;
    const double moment_of_inertia = fibre_length * fibre_length / 12; // rod about its centre, per unit mass
    const double angle = stickiness * torque_magnitude / moment_of_inertia;
    if (angle < 1e-16) {
        return result;
    }
    for (unsigned int i = 0; i < 3; i++) {
        axis[i] /= torque_magnitude;
    }

    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const std::vector<double> axis_cross_direction = cross(axis, old_direction);

    std::vector<double> new_direction(3, 0.0);
    double new_norm = 1e-32; // same guard as BioFVM::normalize
    for (unsigned int i = 0; i < 3; i++) {
        new_direction[i] = old_direction[i] * c + axis_cross_direction[i] * s;
        new_norm += new_direction[i] * new_direction[i];
    }
    new_norm = std::sqrt(new_norm);
    for (unsigned int i = 0; i < 3; i++) {
        new_direction[i] /= new_norm;
    }

    std::vector<double> pivot_to_centre(3, 0.0);
    for (unsigned int i = 0; i < 3; i++) {
        pivot_to_centre[i] = centre[i] - pivot[i];
    }
    const double pivot_to_centre_distance = norm(pivot_to_centre);

    result.rotated = true;
    result.orientation = new_direction;
    result.centre = centre;
    if (pivot_to_centre_distance > 1e-10) {
        // the centre is on the fibre axis, on the side of the impact point or on the other side of the pivot
        const double side = (dot(pivot_to_centre, old_direction) < 0) ? -1.0 : 1.0;
        for (unsigned int i = 0; i < 3; i++) {
            result.centre[i] += side * (-old_direction[i] + new_direction[i]) * pivot_to_centre_distance;
        }
    }
    return result;
}

} // namespace PhysiMeSS_Geometry

#endif
