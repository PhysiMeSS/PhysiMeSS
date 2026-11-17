/*
    Unit tests of the geometry of the cell-fibre interaction of PhysiMeSS
    (addons/PhysiMeSS/PhysiMeSS_geometry.h).

    The expected values are derived from the mathematics alone: closed forms
    for simple configurations, and properties which must hold for any
    configuration (rigid rotation, invariance under a change of frame, ...).
    The rotations used to build the expected values are computed here with
    the full Rodrigues formula, independently of the code under test.

    No test framework: "make physimess-tests" from the PhysiCell root, or
    "make run" in this directory. The exit code is the number of failed tests.
*/

#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "../../../addons/PhysiMeSS/PhysiMeSS_geometry.h"

using namespace PhysiMeSS_Geometry;

typedef std::vector<double> Vec;

// ---------------------------------------------------------------- harness

static int checks_failed_in_test = 0;

#define CHECK(condition) \
    do { \
        if (!(condition)) { \
            checks_failed_in_test++; \
            std::printf("    line %d: CHECK(%s) failed\n", __LINE__, #condition); \
        } \
    } while (0)

#define CHECK_NEAR(value, expected, tolerance) \
    do { \
        const double v_ = (value), e_ = (expected); \
        if (!(std::fabs(v_ - e_) <= (tolerance))) { \
            checks_failed_in_test++; \
            std::printf("    line %d: %s = %.17g, expected %.17g\n", __LINE__, #value, v_, e_); \
        } \
    } while (0)

#define CHECK_VEC_NEAR(value, expected, tolerance) \
    do { \
        const Vec v_ = (value), e_ = (expected); \
        if (!(norm(sub(v_, e_)) <= (tolerance))) { \
            checks_failed_in_test++; \
            std::printf("    line %d: %s = (%.17g, %.17g, %.17g), expected (%.17g, %.17g, %.17g)\n", \
                        __LINE__, #value, v_[0], v_[1], v_[2], e_[0], e_[1], e_[2]); \
        } \
    } while (0)

static const double TOL = 1e-12;
static const double PI = 3.14159265358979323846;

// ---------------------------------------------------------------- helpers

static Vec vec( double x, double y, double z )
{
    Vec v(3);
    v[0] = x; v[1] = y; v[2] = z;
    return v;
}

static Vec add( const Vec& a, const Vec& b ) { return vec(a[0]+b[0], a[1]+b[1], a[2]+b[2]); }
static Vec sub( const Vec& a, const Vec& b ) { return vec(a[0]-b[0], a[1]-b[1], a[2]-b[2]); }
static Vec scale( double s, const Vec& a ) { return vec(s*a[0], s*a[1], s*a[2]); }
static Vec unit( const Vec& a ) { return scale(1.0 / norm(a), a); }

// rotation of v by `angle` about the unit vector k (full Rodrigues formula)
static Vec rotate( const Vec& v, const Vec& k, double angle )
{
    const double c = std::cos(angle), s = std::sin(angle);
    return add(add(scale(c, v), scale(s, cross(k, v))), scale((1 - c) * dot(k, v), k));
}

static std::mt19937 generator(20240613);

static double uniform( double low, double high )
{
    return low + (high - low) * (generator() / 4294967296.0);
}

static Vec random_vector( double size )
{
    return vec(uniform(-size, size), uniform(-size, size), uniform(-size, size));
}

static Vec random_unit_vector()
{
    Vec v = random_vector(1.0);
    while (norm(v) < 0.1 || norm(v) > 1.0) {
        v = random_vector(1.0);
    }
    return unit(v);
}

// a unit vector perpendicular to u
static Vec random_perpendicular( const Vec& u )
{
    Vec v = cross(u, random_unit_vector());
    while (norm(v) < 0.1) {
        v = cross(u, random_unit_vector());
    }
    return unit(v);
}

// ------------------------------------------------- orientation of a fibre

static void orientation_from_angles_closed_forms()
{
    CHECK_VEC_NEAR(orientation_from_angles(0, PI / 2), vec(1, 0, 0), TOL);
    CHECK_VEC_NEAR(orientation_from_angles(PI / 2, PI / 2), vec(0, 1, 0), TOL);
    CHECK_VEC_NEAR(orientation_from_angles(PI, PI / 2), vec(-1, 0, 0), TOL);
    // along z whatever theta
    CHECK_VEC_NEAR(orientation_from_angles(0, 0), vec(0, 0, 1), TOL);
    CHECK_VEC_NEAR(orientation_from_angles(1.234, 0), vec(0, 0, 1), TOL);
    CHECK_VEC_NEAR(orientation_from_angles(1.234, PI), vec(0, 0, -1), TOL);
    // 45 degrees from the z axis, above the diagonal of the xy plane
    CHECK_VEC_NEAR(orientation_from_angles(PI / 4, PI / 4), vec(0.5, 0.5, std::sqrt(0.5)), TOL);
}

static void orientation_from_angles_properties()
{
    for (int n = 0; n < 500; n++) {
        const double theta = uniform(-PI, PI), phi = uniform(0.01, PI - 0.01);
        const Vec orientation = orientation_from_angles(theta, phi);
        CHECK_NEAR(norm(orientation), 1.0, TOL);
        // phi is the angle from the z axis, theta the angle of the projection on the xy plane
        CHECK_NEAR(std::acos(orientation[2]), phi, 1e-9);
        CHECK_NEAR(std::atan2(orientation[1], orientation[0]), theta, 1e-9);
        // with phi = pi/2: the orientation of 2D simulations
        const Vec planar = orientation_from_angles(theta, PI / 2);
        CHECK_VEC_NEAR(planar, vec(std::cos(theta), std::sin(theta), 0), TOL);
    }
}

// ------------------------------------------------- nearest point on a fibre

static void nearest_point_along_the_fibre()
{
    // fibre from (-10,0,0) to (10,0,0)
    const Vec centre = vec(0, 0, 0), orientation = vec(1, 0, 0);
    CHECK_VEC_NEAR(displacement_from_fibre(vec(3, 4, 0), centre, orientation, 10), vec(0, 4, 0), TOL);
    CHECK_VEC_NEAR(displacement_from_fibre(vec(-7, 0, -2), centre, orientation, 10), vec(0, 0, -2), TOL);
    CHECK_VEC_NEAR(displacement_from_fibre(vec(0, 1, 1), centre, orientation, 10), vec(0, 1, 1), TOL);
}

static void nearest_point_beyond_the_ends()
{
    const Vec centre = vec(0, 0, 0), orientation = vec(1, 0, 0);
    CHECK_VEC_NEAR(displacement_from_fibre(vec(13, 4, 0), centre, orientation, 10), vec(3, 4, 0), TOL);
    CHECK_VEC_NEAR(displacement_from_fibre(vec(-13, 4, 0), centre, orientation, 10), vec(-3, 4, 0), TOL);
    CHECK_VEC_NEAR(displacement_from_fibre(vec(15, 0, 0), centre, orientation, 10), vec(5, 0, 0), TOL);
    CHECK_VEC_NEAR(displacement_from_fibre(vec(-10, 0, 6), centre, orientation, 10), vec(0, 0, 6), TOL);
}

static void nearest_point_of_a_point_on_the_fibre()
{
    const Vec centre = vec(1, 2, 3), orientation = unit(vec(1, 1, 1));
    const double half_length = 5;
    for (int i = -4; i <= 4; i++) {
        const Vec point = add(centre, scale(i * half_length / 4.0, orientation));
        CHECK_NEAR(norm(displacement_from_fibre(point, centre, orientation, half_length)), 0.0, TOL);
    }
}

static void nearest_point_oblique_fibre()
{
    const Vec centre = vec(1, 2, 3), orientation = unit(vec(1, 1, 1));
    const Vec normal = unit(vec(1, -1, 0)); // perpendicular to the fibre
    const double half_length = 5;

    // above a point of the fibre
    Vec point = add(add(centre, scale(2, orientation)), scale(3, normal));
    CHECK_VEC_NEAR(displacement_from_fibre(point, centre, orientation, half_length), scale(3, normal), TOL);

    // beyond the end: the nearest point is the end
    point = add(add(centre, scale(8, orientation)), scale(4, normal));
    CHECK_VEC_NEAR(displacement_from_fibre(point, centre, orientation, half_length),
                   add(scale(3, orientation), scale(4, normal)), TOL);
    CHECK_NEAR(norm(displacement_from_fibre(point, centre, orientation, half_length)), 5.0, TOL);
}

static void nearest_point_does_not_depend_on_the_sign_of_the_orientation()
{
    for (int n = 0; n < 200; n++) {
        const Vec centre = random_vector(50), orientation = random_unit_vector(), point = random_vector(80);
        const double half_length = uniform(1, 40);
        CHECK_VEC_NEAR(displacement_from_fibre(point, centre, orientation, half_length),
                       displacement_from_fibre(point, centre, scale(-1, orientation), half_length), 1e-10);
    }
}

static void nearest_point_is_on_the_fibre_and_is_the_nearest()
{
    for (int n = 0; n < 500; n++) {
        const Vec centre = random_vector(50), orientation = random_unit_vector(), point = random_vector(80);
        const double half_length = uniform(1, 40);
        const Vec displacement = displacement_from_fibre(point, centre, orientation, half_length);
        const Vec nearest = sub(point, displacement);

        // on the fibre axis, between the two ends
        const Vec from_centre = sub(nearest, centre);
        const double abscissa = dot(from_centre, orientation);
        CHECK_NEAR(norm(cross(from_centre, orientation)), 0.0, 1e-10);
        CHECK(std::fabs(abscissa) <= half_length + 1e-10);

        // perpendicular to the fibre unless the nearest point is an end
        if (std::fabs(abscissa) < half_length - 1e-8) {
            CHECK_NEAR(dot(displacement, orientation), 0.0, 1e-10);
        }

        // no point of the fibre is nearer
        for (int i = -100; i <= 100; i++) {
            const Vec other = add(centre, scale(i * half_length / 100.0, orientation));
            CHECK(norm(displacement) <= norm(sub(point, other)) + 1e-10);
        }
    }
}

static void nearest_point_does_not_depend_on_the_frame()
{
    for (int n = 0; n < 200; n++) {
        const Vec centre = random_vector(50), orientation = random_unit_vector(), point = random_vector(80);
        const double half_length = uniform(1, 40);
        const Vec k = random_unit_vector(), shift = random_vector(100);
        const double angle = uniform(-PI, PI);

        const Vec reference = displacement_from_fibre(point, centre, orientation, half_length);
        const Vec moved = displacement_from_fibre(add(rotate(point, k, angle), shift),
                                                  add(rotate(centre, k, angle), shift),
                                                  rotate(orientation, k, angle), half_length);
        CHECK_VEC_NEAR(moved, rotate(reference, k, angle), 1e-9);
    }
}

static void nearest_point_of_a_fibre_of_zero_length()
{
    // a fibre reduced to a point: the displacement is point - centre
    const Vec displacement = displacement_from_fibre(vec(4, 5, 6), vec(1, 1, 1), vec(1, 0, 0), 0.0);
    CHECK_VEC_NEAR(displacement, vec(3, 4, 5), TOL);
}

// ------------------------------------- angular velocity given by a cell
// angular_velocity_about_pivot(pivot, centre, impact_point, velocity, half_length)

static void angular_velocity_of_a_free_fibre()
{
    // fibre from (-10,0,0) to (10,0,0), free: the pivot is the centre, L^2/12 = 100/3
    const Vec centre = vec(0, 0, 0);
    // dragged towards +y at x = +8: anticlockwise seen from +z
    CHECK_VEC_NEAR(angular_velocity_about_pivot(centre, centre, vec(8, 0, 0), vec(0, 0.5, 0), 10), vec(0, 0, 0.12), TOL);
    // dragged towards +y at x = -8: clockwise seen from +z
    CHECK_VEC_NEAR(angular_velocity_about_pivot(centre, centre, vec(-8, 0, 0), vec(0, 0.5, 0), 10), vec(0, 0, -0.12), TOL);
    // dragged towards +z at x = +8: out of the plane, about -y
    CHECK_VEC_NEAR(angular_velocity_about_pivot(centre, centre, vec(8, 0, 0), vec(0, 0, 0.5), 10), vec(0, -0.12, 0), TOL);
    // fibre along z dragged towards +x at z = -8: about -y as well
    CHECK_VEC_NEAR(angular_velocity_about_pivot(centre, centre, vec(0, 0, -8), vec(0.5, 0, 0), 10), vec(0, -0.12, 0), TOL);
    // oblique velocity (1, 2, 2) at x = +6, fibre of length 24 (L^2/12 = 48): (6,0,0) x (1,2,2) = (0,-12,12)
    CHECK_VEC_NEAR(angular_velocity_about_pivot(centre, centre, vec(6, 0, 0), vec(1, 2, 2), 12), vec(0, -0.25, 0.25), TOL);
}

static void angular_velocity_of_a_crosslinked_fibre()
{
    // fibre from (0,0,0) to (20,0,0) crosslinked at its end: L^2/12 + d^2 = 100/3 + 100 = L^2/3
    CHECK_VEC_NEAR(angular_velocity_about_pivot(vec(0, 0, 0), vec(10, 0, 0), vec(15, 0, 0), vec(0, 0.4, 0), 10),
                   vec(0, 0, 15 * 0.4 / (400.0 / 3)), TOL);
    // fibre from (-10,0,0) to (10,0,0) crosslinked at x = -5: 100/3 + 25
    CHECK_VEC_NEAR(angular_velocity_about_pivot(vec(-5, 0, 0), vec(0, 0, 0), vec(5, 0, 0), vec(0, 0, 0.4), 10),
                   vec(0, -10 * 0.4 / (100.0 / 3 + 25), 0), TOL);
    // ... dragged on the other side of the crosslink, at x = -9: the rotation is reversed
    CHECK_VEC_NEAR(angular_velocity_about_pivot(vec(-5, 0, 0), vec(0, 0, 0), vec(-9, 0, 0), vec(0, 0, 0.4), 10),
                   vec(0, 4 * 0.4 / (100.0 / 3 + 25), 0), TOL);

    // with the same lever arm, a fibre is 4 times harder to rotate about its end than about its centre
    const Vec about_centre = angular_velocity_about_pivot(vec(0, 0, 0), vec(0, 0, 0), vec(5, 0, 0), vec(0, 1, 0), 10);
    const Vec about_end = angular_velocity_about_pivot(vec(-10, 0, 0), vec(0, 0, 0), vec(-5, 0, 0), vec(0, 1, 0), 10);
    CHECK_NEAR(norm(about_centre) / norm(about_end), 4.0, 1e-10);
}

static void angular_velocity_balances_the_drag_along_the_fibre()
{
    // Independent derivation of L^2/12 + d^2. Overdamped rod with a drag coefficient of 1
    // per unit length: a force F translates it at the velocity v = F / L. The same force
    // applied at the lever arm r exerts the torque r x F about the pivot, balanced by the
    // drag of all the points x of the rod, which move at omega x (x - pivot):
    //     sum over the rod of  (x - pivot) x (omega x (x - pivot)) ds  =  r x (L v)
    for (int n = 0; n < 200; n++) {
        const Vec centre = random_vector(50), direction = random_unit_vector();
        const double half_length = uniform(5, 40), length = 2 * half_length;
        const Vec pivot = add(centre, scale(uniform(-1, 1) * half_length, direction));
        const Vec impact = add(centre, scale(uniform(-1, 1) * half_length, direction));
        const Vec velocity = random_vector(1.0);

        const Vec omega = angular_velocity_about_pivot(pivot, centre, impact, velocity, half_length);

        const int pieces = 2000;
        const double ds = length / pieces;
        Vec drag_torque = vec(0, 0, 0);
        for (int i = 0; i < pieces; i++) {
            const Vec x = sub(add(centre, scale(-half_length + (i + 0.5) * ds, direction)), pivot);
            drag_torque = add(drag_torque, scale(ds, cross(x, cross(omega, x))));
        }
        const Vec applied_torque = cross(sub(impact, pivot), scale(length, velocity));
        CHECK_VEC_NEAR(drag_torque, applied_torque, 1e-5 * (1 + norm(applied_torque)));
    }
}

static void no_angular_velocity_without_torque()
{
    const Vec centre = vec(1, 2, 3), direction = unit(vec(1, 2, -2));
    const Vec impact = add(centre, scale(4, direction));
    const Vec crosslink = add(centre, scale(-6, direction));
    // cell at the pivot
    CHECK_NEAR(norm(angular_velocity_about_pivot(centre, centre, centre, vec(0, 1, 0), 10)), 0.0, TOL);
    CHECK_NEAR(norm(angular_velocity_about_pivot(crosslink, centre, crosslink, vec(0, 1, 0), 10)), 0.0, TOL);
    // cell moving along the fibre
    CHECK_NEAR(norm(angular_velocity_about_pivot(centre, centre, impact, scale(2, direction), 10)), 0.0, TOL);
    CHECK_NEAR(norm(angular_velocity_about_pivot(crosslink, centre, impact, scale(-2, direction), 10)), 0.0, TOL);
    // cell at rest
    CHECK_NEAR(norm(angular_velocity_about_pivot(centre, centre, impact, vec(0, 0, 0), 10)), 0.0, TOL);
}

static void angular_velocities_of_several_cells_add_up()
{
    const Vec centre = vec(0, 0, 0);
    // two cells dragging the two ends in the same direction: no rotation
    Vec total = add(angular_velocity_about_pivot(centre, centre, vec(8, 0, 0), vec(0, 0.5, 0), 10),
                    angular_velocity_about_pivot(centre, centre, vec(-8, 0, 0), vec(0, 0.5, 0), 10));
    CHECK_NEAR(norm(total), 0.0, TOL);
    // two cells dragging the two ends in opposite directions: twice the rotation
    total = add(angular_velocity_about_pivot(centre, centre, vec(8, 0, 0), vec(0, 0.5, 0), 10),
                angular_velocity_about_pivot(centre, centre, vec(-8, 0, 0), vec(0, -0.5, 0), 10));
    CHECK_VEC_NEAR(total, vec(0, 0, 0.24), TOL);

    // linear in the velocity and in the lever arm
    for (int n = 0; n < 200; n++) {
        const Vec direction = random_unit_vector(), pivot = random_vector(50);
        const double half_length = uniform(5, 40);
        const Vec fibre_centre = add(pivot, scale(uniform(-1, 1) * half_length, direction));
        const Vec impact = add(pivot, scale(uniform(-10, 10), direction));
        const Vec v1 = random_vector(1.0), v2 = random_vector(1.0);
        const double a = uniform(-2, 2);

        const Vec w1 = angular_velocity_about_pivot(pivot, fibre_centre, impact, v1, half_length);
        const Vec w2 = angular_velocity_about_pivot(pivot, fibre_centre, impact, v2, half_length);
        CHECK_VEC_NEAR(angular_velocity_about_pivot(pivot, fibre_centre, impact, add(v1, scale(a, v2)), half_length),
                       add(w1, scale(a, w2)), 1e-10);
        const Vec twice_as_far = add(pivot, scale(2, sub(impact, pivot)));
        CHECK_VEC_NEAR(angular_velocity_about_pivot(pivot, fibre_centre, twice_as_far, v1, half_length), scale(2, w1), 1e-10);
    }
}

static void angular_velocity_properties_in_3D()
{
    for (int n = 0; n < 1000; n++) {
        const Vec centre = random_vector(50), direction = random_unit_vector();
        const double half_length = uniform(5, 40);
        const Vec pivot = add(centre, scale(uniform(-1, 1) * half_length, direction));
        const Vec impact = add(centre, scale(uniform(-1, 1) * half_length, direction));
        const Vec lever = sub(impact, pivot);
        const Vec velocity = random_vector(1.0);
        const Vec omega = angular_velocity_about_pivot(pivot, centre, impact, velocity, half_length);

        // perpendicular to the fibre and to the velocity of the cell
        CHECK_NEAR(dot(omega, direction), 0.0, 1e-10);
        CHECK_NEAR(dot(omega, velocity), 0.0, 1e-10);
        // right hand rule: the point of the fibre touched by the cell follows the cell
        CHECK(dot(cross(omega, lever), velocity) >= 0);
        // the velocity of the cell along the fibre has no effect
        CHECK_VEC_NEAR(angular_velocity_about_pivot(pivot, centre, impact, add(velocity, scale(uniform(-5, 5), direction)), half_length),
                       omega, 1e-10);

        // same result in another frame
        const Vec k = random_unit_vector(), shift = random_vector(100);
        const double angle = uniform(-PI, PI);
        CHECK_VEC_NEAR(angular_velocity_about_pivot(add(rotate(pivot, k, angle), shift), add(rotate(centre, k, angle), shift),
                                                    add(rotate(impact, k, angle), shift), rotate(velocity, k, angle), half_length),
                       rotate(omega, k, angle), 1e-10);
    }
}

// ------------------------------------------------- rotation of a fibre
// rotate_fibre_about_pivot(pivot, centre, orientation, rotation)

static void rotation_closed_forms()
{
    const Vec origin = vec(0, 0, 0);
    const double angle = 0.3;
    // fibre along x rotated about +z: anticlockwise in the xy plane
    FibreRotation rotation = rotate_fibre_about_pivot(origin, origin, vec(1, 0, 0), vec(0, 0, angle));
    CHECK(rotation.rotated);
    CHECK_VEC_NEAR(rotation.orientation, vec(std::cos(angle), std::sin(angle), 0), TOL);
    CHECK_VEC_NEAR(rotation.centre, origin, TOL);

    // about -y: the positive end rises towards +z
    rotation = rotate_fibre_about_pivot(origin, origin, vec(1, 0, 0), vec(0, -angle, 0));
    CHECK_VEC_NEAR(rotation.orientation, vec(std::cos(angle), 0, std::sin(angle)), TOL);

    // a third of a turn about (1,1,1) maps x to y
    rotation = rotate_fibre_about_pivot(origin, origin, vec(1, 0, 0), scale(2 * PI / 3, unit(vec(1, 1, 1))));
    CHECK_VEC_NEAR(rotation.orientation, vec(0, 1, 0), TOL);

    // fibre from (0,0,0) to (20,0,0) hinged at the origin, rotated about +z then about -y
    rotation = rotate_fibre_about_pivot(origin, vec(10, 0, 0), vec(1, 0, 0), vec(0, 0, angle));
    CHECK_VEC_NEAR(rotation.orientation, vec(std::cos(angle), std::sin(angle), 0), TOL);
    CHECK_VEC_NEAR(rotation.centre, vec(10 * std::cos(angle), 10 * std::sin(angle), 0), TOL);
    rotation = rotate_fibre_about_pivot(origin, vec(10, 0, 0), vec(1, 0, 0), vec(0, -angle, 0));
    CHECK_VEC_NEAR(rotation.centre, vec(10 * std::cos(angle), 0, 10 * std::sin(angle)), TOL);
}

static void rotation_keeps_the_sign_of_the_orientation()
{
    // the same fibre described by the opposite orientation gives the opposite orientation
    const Vec origin = vec(0, 0, 0);
    const FibreRotation rotation = rotate_fibre_about_pivot(origin, origin, vec(-1, 0, 0), vec(0, 0, 0.3));
    CHECK_VEC_NEAR(rotation.orientation, vec(-std::cos(0.3), -std::sin(0.3), 0), TOL);

    for (int n = 0; n < 200; n++) {
        const Vec pivot = random_vector(50), centre = random_vector(50), orientation = random_unit_vector();
        const Vec turn = scale(uniform(0.001, 1.5), random_unit_vector());
        const FibreRotation a = rotate_fibre_about_pivot(pivot, centre, orientation, turn);
        const FibreRotation b = rotate_fibre_about_pivot(pivot, centre, scale(-1, orientation), turn);
        CHECK_VEC_NEAR(b.orientation, scale(-1, a.orientation), 1e-10);
        CHECK_VEC_NEAR(b.centre, a.centre, 1e-9);
        CHECK(dot(a.orientation, orientation) > 0); // less than a quarter turn
    }
}

static void rotation_is_rigid()
{
    for (int n = 0; n < 1000; n++) {
        const Vec centre = random_vector(50), orientation = random_unit_vector();
        const double half_length = uniform(5, 40);
        // pivot anywhere on the fibre, on either side of the centre
        const double abscissa = uniform(-1, 1) * half_length;
        const Vec pivot = add(centre, scale(abscissa, orientation));
        // any rotation vector, not necessarily perpendicular to the fibre
        const Vec axis = random_unit_vector();
        const double angle = uniform(0.001, PI);

        const FibreRotation rotation = rotate_fibre_about_pivot(pivot, centre, orientation, scale(angle, axis));
        CHECK(rotation.rotated);
        CHECK_NEAR(norm(rotation.orientation), 1.0, TOL);
        CHECK_VEC_NEAR(rotation.orientation, rotate(orientation, axis, angle), 1e-10);
        CHECK_VEC_NEAR(rotation.centre, add(pivot, rotate(sub(centre, pivot), axis, angle)), 1e-9);

        // the pivot stays at the same place on the fibre
        const Vec centre_to_pivot = sub(pivot, rotation.centre);
        CHECK_NEAR(dot(centre_to_pivot, rotation.orientation), abscissa, 1e-9);
        CHECK_NEAR(norm(cross(centre_to_pivot, rotation.orientation)), 0.0, 1e-9);
    }
}

static void rotation_about_the_fibre_axis_changes_nothing()
{
    for (int n = 0; n < 100; n++) {
        const Vec centre = random_vector(50), orientation = random_unit_vector();
        const Vec pivot = add(centre, scale(uniform(-20, 20), orientation));
        const FibreRotation rotation = rotate_fibre_about_pivot(pivot, centre, orientation, scale(uniform(-3, 3), orientation));
        if (rotation.rotated) {
            CHECK_VEC_NEAR(rotation.orientation, orientation, 1e-10);
            CHECK_VEC_NEAR(rotation.centre, centre, 1e-9);
        }
    }
}

static void no_rotation_for_a_zero_rotation_vector()
{
    CHECK(!rotate_fibre_about_pivot(vec(0, 0, 0), vec(1, 2, 3), vec(1, 0, 0), vec(0, 0, 0)).rotated);
    CHECK(!rotate_fibre_about_pivot(vec(0, 0, 0), vec(1, 2, 3), vec(1, 0, 0), vec(0, 1e-18, 0)).rotated);
}

static void rotation_in_small_steps_equals_one_rotation()
{
    // an angular velocity applied during dt in N time steps of dt/N gives the same fibre
    for (int n = 0; n < 100; n++) {
        const Vec pivot = random_vector(50), orientation = random_unit_vector();
        const Vec centre = add(pivot, scale(uniform(-20, 20), orientation));
        const Vec turn = scale(uniform(0.01, 2), random_unit_vector());

        const FibreRotation once = rotate_fibre_about_pivot(pivot, centre, orientation, turn);

        const int steps = 100;
        Vec stepped_centre = centre, stepped_orientation = orientation;
        for (int i = 0; i < steps; i++) {
            const FibreRotation step = rotate_fibre_about_pivot(pivot, stepped_centre, stepped_orientation, scale(1.0 / steps, turn));
            CHECK(step.rotated);
            stepped_centre = step.centre;
            stepped_orientation = step.orientation;
        }
        CHECK_VEC_NEAR(stepped_orientation, once.orientation, 1e-9);
        CHECK_VEC_NEAR(stepped_centre, once.centre, 1e-8);
    }
}

static void rotation_in_the_xy_plane_matches_the_2D_rotation()
{
    for (int n = 0; n < 500; n++) {
        const double phi = uniform(-PI, PI);
        const Vec orientation = vec(std::cos(phi), std::sin(phi), 0);
        const Vec pivot = vec(uniform(-50, 50), uniform(-50, 50), 0);
        const Vec centre = add(pivot, scale(uniform(-20, 20), orientation));
        const double angle = uniform(-1.5, 1.5); // signed angle about +z

        const FibreRotation rotation = rotate_fibre_about_pivot(pivot, centre, orientation, vec(0, 0, angle));
        CHECK(rotation.rotated);
        if (!rotation.rotated) {
            continue;
        }
        // 2D rotation matrix [c -s; s c]
        const double c = std::cos(angle), s = std::sin(angle);
        const Vec arm = sub(centre, pivot);
        CHECK_VEC_NEAR(rotation.orientation, vec(orientation[0] * c - orientation[1] * s, orientation[0] * s + orientation[1] * c, 0), 1e-10);
        CHECK_VEC_NEAR(rotation.centre, add(pivot, vec(arm[0] * c - arm[1] * s, arm[0] * s + arm[1] * c, 0)), 1e-9);
        // stays exactly in the plane
        CHECK(rotation.orientation[2] == 0.0);
        CHECK(rotation.centre[2] == 0.0);
    }
}

static void rotation_does_not_depend_on_the_frame()
{
    for (int n = 0; n < 500; n++) {
        const Vec pivot = random_vector(50), centre = random_vector(50), orientation = random_unit_vector();
        const Vec turn = scale(uniform(0.001, 2), random_unit_vector());
        const Vec k = random_unit_vector(), shift = random_vector(100);
        const double angle = uniform(-PI, PI);

        const FibreRotation reference = rotate_fibre_about_pivot(pivot, centre, orientation, turn);
        const FibreRotation moved = rotate_fibre_about_pivot(add(rotate(pivot, k, angle), shift), add(rotate(centre, k, angle), shift),
                                                             rotate(orientation, k, angle), rotate(turn, k, angle));
        CHECK(reference.rotated && moved.rotated);
        CHECK_VEC_NEAR(moved.orientation, rotate(reference.orientation, k, angle), 1e-9);
        CHECK_VEC_NEAR(moved.centre, add(rotate(reference.centre, k, angle), shift), 1e-9);
    }
}

// ------------------------------------- a cell dragging a fibre during one time step
// rotation vector = dt * angular velocity of the traction, fibre_sticky * motility_vector

static FibreRotation cell_rotates_fibre( const Vec& pivot, const Vec& centre, const Vec& orientation, const Vec& impact,
                                         const Vec& velocity, double half_length, double stickiness, double dt )
{
    const Vec omega = angular_velocity_about_pivot(pivot, centre, impact, scale(stickiness, velocity), half_length);
    return rotate_fibre_about_pivot(pivot, centre, orientation, scale(dt, omega));
}

static void free_fibre_dragged_by_a_cell()
{
    const Vec centre = vec(0, 0, 0), orientation = vec(1, 0, 0);
    const double half_length = 10, stickiness = 2.0, dt = 0.1;
    const double angle = stickiness * dt * 0.12;

    // towards +y at the positive end: anticlockwise
    FibreRotation rotation = cell_rotates_fibre(centre, centre, orientation, vec(8, 0, 0), vec(0, 0.5, 0), half_length, stickiness, dt);
    CHECK(rotation.rotated);
    CHECK_VEC_NEAR(rotation.orientation, vec(std::cos(angle), std::sin(angle), 0), TOL);
    CHECK_VEC_NEAR(rotation.centre, centre, TOL);

    // towards +y at the negative end: clockwise, and the orientation is not flipped
    rotation = cell_rotates_fibre(centre, centre, orientation, vec(-8, 0, 0), vec(0, 0.5, 0), half_length, stickiness, dt);
    CHECK_VEC_NEAR(rotation.orientation, vec(std::cos(angle), -std::sin(angle), 0), TOL);

    // towards +z at the positive end: out of the plane
    rotation = cell_rotates_fibre(centre, centre, orientation, vec(8, 0, 0), vec(0, 0, 0.5), half_length, stickiness, dt);
    CHECK_VEC_NEAR(rotation.orientation, vec(std::cos(angle), 0, std::sin(angle)), TOL);

    // half the time step, half the angle; no stickiness, no rotation
    rotation = cell_rotates_fibre(centre, centre, orientation, vec(8, 0, 0), vec(0, 0.5, 0), half_length, stickiness, dt / 2);
    CHECK_VEC_NEAR(rotation.orientation, vec(std::cos(angle / 2), std::sin(angle / 2), 0), TOL);
    CHECK(!cell_rotates_fibre(centre, centre, orientation, vec(8, 0, 0), vec(0, 0.5, 0), half_length, 0.0, dt).rotated);
}

static void crosslinked_fibre_dragged_by_a_cell()
{
    const double half_length = 10, stickiness = 2.0, dt = 0.1;

    // fibre from (0,0,0) to (20,0,0) hinged at the origin, dragged towards +y at x = 15
    double angle = stickiness * dt * 15 * 0.4 / (400.0 / 3);
    FibreRotation rotation = cell_rotates_fibre(vec(0, 0, 0), vec(10, 0, 0), vec(1, 0, 0), vec(15, 0, 0), vec(0, 0.4, 0),
                                                half_length, stickiness, dt);
    CHECK(rotation.rotated);
    CHECK_VEC_NEAR(rotation.orientation, vec(std::cos(angle), std::sin(angle), 0), TOL);
    CHECK_VEC_NEAR(rotation.centre, vec(10 * std::cos(angle), 10 * std::sin(angle), 0), TOL);

    // same fibre dragged towards +z: it leaves the plane, hinged at the origin
    rotation = cell_rotates_fibre(vec(0, 0, 0), vec(10, 0, 0), vec(1, 0, 0), vec(15, 0, 0), vec(0, 0, 0.4), half_length, stickiness, dt);
    CHECK_VEC_NEAR(rotation.orientation, vec(std::cos(angle), 0, std::sin(angle)), TOL);
    CHECK_VEC_NEAR(rotation.centre, vec(10 * std::cos(angle), 0, 10 * std::sin(angle)), TOL);

    // fibre from (-10,0,0) to (10,0,0) crosslinked at x = -4, dragged towards +y at x = -9, on the
    // side of the crosslink opposite to the centre: clockwise, the centre (at 4 from the crosslink) goes down
    angle = stickiness * dt * 5 * 0.4 / (100.0 / 3 + 16);
    rotation = cell_rotates_fibre(vec(-4, 0, 0), vec(0, 0, 0), vec(1, 0, 0), vec(-9, 0, 0), vec(0, 0.4, 0), half_length, stickiness, dt);
    CHECK(rotation.rotated);
    CHECK_VEC_NEAR(rotation.orientation, vec(std::cos(angle), -std::sin(angle), 0), TOL);
    CHECK_VEC_NEAR(rotation.centre, vec(-4 + 4 * std::cos(angle), -4 * std::sin(angle), 0), TOL);
}

static void the_point_touched_by_the_cell_follows_the_cell()
{
    for (int n = 0; n < 1000; n++) {
        const Vec centre = random_vector(50), orientation = random_unit_vector();
        const double half_length = uniform(5, 40);
        // free fibre, or crosslinked anywhere along its length
        const Vec pivot = (n % 2 == 0) ? centre : add(centre, scale(uniform(-1, 1) * half_length, orientation));
        const double abscissa = uniform(-1, 1) * half_length;
        const Vec impact = add(centre, scale(abscissa, orientation));
        if (norm(sub(impact, pivot)) < 0.1) {
            continue;
        }
        const Vec velocity = add(scale(uniform(0.05, 0.5), random_perpendicular(orientation)), scale(uniform(-1, 1), orientation));

        const FibreRotation rotation = cell_rotates_fibre(pivot, centre, orientation, impact, velocity, half_length, 1.0, 0.1);
        CHECK(rotation.rotated);
        if (!rotation.rotated) {
            continue;
        }
        // the same material point of the fibre, after the rotation
        const Vec new_impact = add(rotation.centre, scale(abscissa, rotation.orientation));
        CHECK(dot(sub(new_impact, impact), velocity) > 0);
        // it stays at the same distance from the pivot
        CHECK_NEAR(norm(sub(new_impact, pivot)), norm(sub(impact, pivot)), 1e-9);
    }
}

// ------------------------------------- contact between a cell and a fibre
// contact_velocity(displacement, contact_distance, repulsion)

static void contact_velocity_closed_forms()
{
    // cell centre at 3 above the fibre axis, contact up to 6 (cell radius + fibre radius), repulsion 0.4:
    // the fibre is pushed away from the cell at 0.4 * (1 - 3/6)^2 = 0.1
    CHECK_VEC_NEAR(contact_velocity(vec(0, 3, 0), 6, 0.4), vec(0, -0.1, 0), TOL);
    CHECK_VEC_NEAR(contact_velocity(vec(0, 0, -3), 6, 0.4), vec(0, 0, 0.1), TOL);
    // oblique: displacement (2, -1, 2) of norm 3
    CHECK_VEC_NEAR(contact_velocity(vec(2, -1, 2), 6, 0.4), scale(-0.1 / 3, vec(2, -1, 2)), TOL);
    // no contact: at the contact distance or beyond
    CHECK_NEAR(norm(contact_velocity(vec(0, 6, 0), 6, 0.4)), 0.0, TOL);
    CHECK_NEAR(norm(contact_velocity(vec(0, 7, 0), 6, 0.4)), 0.0, TOL);
    // no repulsion
    CHECK_NEAR(norm(contact_velocity(vec(0, 3, 0), 6, 0.0)), 0.0, TOL);
    // cell centre on the fibre axis: no direction to push to, but no NaN
    const Vec on_axis = contact_velocity(vec(0, 0, 0), 6, 0.4);
    CHECK(on_axis[0] == on_axis[0] && on_axis[1] == on_axis[1] && on_axis[2] == on_axis[2]);
    CHECK_NEAR(norm(on_axis), 0.0, TOL);
}

static void contact_velocity_properties()
{
    for (int n = 0; n < 500; n++) {
        const double contact_distance = uniform(5, 15), repulsion = uniform(0.1, 20);
        const Vec direction = random_unit_vector();
        const double d = uniform(0.01, 1.2) * contact_distance;
        const Vec velocity = contact_velocity(scale(d, direction), contact_distance, repulsion);

        if (d >= contact_distance) {
            CHECK_NEAR(norm(velocity), 0.0, TOL);
            continue;
        }
        // along the line of the centres, away from the cell
        CHECK_NEAR(norm(cross(velocity, direction)), 0.0, 1e-10);
        CHECK(dot(velocity, direction) <= 0);
        // PhysiCell repulsion: repulsion * (1 - d/R)^2
        const double overlap = 1 - d / contact_distance;
        CHECK_NEAR(norm(velocity), repulsion * overlap * overlap, 1e-10);
        // the deeper the contact, the stronger
        const Vec deeper = contact_velocity(scale(d / 2, direction), contact_distance, repulsion);
        CHECK(norm(deeper) > norm(velocity));
    }
}

// ------------------------------------- torque of the contact and of the traction
// as in PhysiMeSS_Fibre: the angular velocity is that of the velocity
//     contact_velocity + fibre_sticky * motility_vector
// applied at the point of the fibre nearest to the cell

static Vec cell_angular_velocity( const Vec& pivot, const Vec& centre, const Vec& orientation, double half_length,
                                  const Vec& cell_position, const Vec& motility, double contact_distance,
                                  double repulsion, double stickiness )
{
    const Vec displacement = displacement_from_fibre(cell_position, centre, orientation, half_length);
    if (norm(displacement) > contact_distance) {
        return vec(0, 0, 0);
    }
    const Vec impact = sub(cell_position, displacement);
    const Vec velocity = add(contact_velocity(displacement, contact_distance, repulsion), scale(stickiness, motility));
    return angular_velocity_about_pivot(pivot, centre, impact, velocity, half_length);
}

static void cell_at_rest_pressing_on_a_fibre()
{
    // fibre from (-10,0,0) to (10,0,0); contact distance 6, repulsion 0.4: a cell at 3 from the axis pushes at 0.1
    const Vec centre = vec(0, 0, 0), orientation = vec(1, 0, 0), at_rest = vec(0, 0, 0);
    const double half_length = 10, R = 6, repulsion = 0.4, drag = 100.0 / 3;

    // above the positive end: the end goes down, clockwise seen from +z
    CHECK_VEC_NEAR(cell_angular_velocity(centre, centre, orientation, half_length, vec(8, 3, 0), at_rest, R, repulsion, 1.0),
                   vec(0, 0, -8 * 0.1 / drag), TOL);
    // above the negative end: anticlockwise
    CHECK_VEC_NEAR(cell_angular_velocity(centre, centre, orientation, half_length, vec(-8, 3, 0), at_rest, R, repulsion, 1.0),
                   vec(0, 0, 8 * 0.1 / drag), TOL);
    // above the centre: pushed, not rotated
    CHECK_NEAR(norm(cell_angular_velocity(centre, centre, orientation, half_length, vec(0, 3, 0), at_rest, R, repulsion, 1.0)), 0.0, TOL);
    // out of the plane, at z = +3 over the positive end: the end goes towards -z, rotation about +y
    CHECK_VEC_NEAR(cell_angular_velocity(centre, centre, orientation, half_length, vec(8, 0, 3), at_rest, R, repulsion, 1.0),
                   vec(0, 8 * 0.1 / drag, 0), TOL);
    // too far to touch
    CHECK_NEAR(norm(cell_angular_velocity(centre, centre, orientation, half_length, vec(8, 7, 0), at_rest, R, repulsion, 1.0)), 0.0, TOL);
    // fibre hinged at its negative end: same torque about the pivot for a longer lever, larger drag (L^2/3)
    CHECK_VEC_NEAR(cell_angular_velocity(vec(-10, 0, 0), centre, orientation, half_length, vec(8, 3, 0), at_rest, R, repulsion, 1.0),
                   vec(0, 0, -18 * 0.1 / (400.0 / 3)), TOL);
}

static void cells_pressing_on_both_sides_of_a_fibre()
{
    const Vec centre = vec(0, 0, 0), orientation = vec(1, 0, 0), at_rest = vec(0, 0, 0);
    const double half_length = 10, R = 6, repulsion = 0.4, drag = 100.0 / 3;
    const Vec above_positive_end = cell_angular_velocity(centre, centre, orientation, half_length, vec(8, 3, 0), at_rest, R, repulsion, 1.0);

    // same side, opposite ends: the fibre is pushed, the torques cancel
    Vec total = add(above_positive_end,
                    cell_angular_velocity(centre, centre, orientation, half_length, vec(-8, 3, 0), at_rest, R, repulsion, 1.0));
    CHECK_NEAR(norm(total), 0.0, TOL);
    // opposite sides, opposite ends: the torques add up
    total = add(above_positive_end,
                cell_angular_velocity(centre, centre, orientation, half_length, vec(-8, -3, 0), at_rest, R, repulsion, 1.0));
    CHECK_VEC_NEAR(total, vec(0, 0, -2 * 8 * 0.1 / drag), TOL);
    // opposite sides, same place: the fibre is squeezed, no rotation
    total = add(above_positive_end,
                cell_angular_velocity(centre, centre, orientation, half_length, vec(8, -3, 0), at_rest, R, repulsion, 1.0));
    CHECK_NEAR(norm(total), 0.0, TOL);
}

static void contact_and_traction_add_up()
{
    const Vec centre = vec(0, 0, 0), orientation = vec(1, 0, 0);
    const double half_length = 10, R = 6, repulsion = 0.4, drag = 100.0 / 3;
    const Vec cell = vec(8, 3, 0); // above the positive end, pushes it down at 0.1

    // cell moving towards the fibre (-y) at 0.5, stickiness 2: contact and traction in the same direction
    CHECK_VEC_NEAR(cell_angular_velocity(centre, centre, orientation, half_length, cell, vec(0, -0.5, 0), R, repulsion, 2.0),
                   vec(0, 0, -8 * (0.1 + 2 * 0.5) / drag), TOL);
    // cell moving away from the fibre (+y): the traction pulls the fibre, against the contact
    CHECK_VEC_NEAR(cell_angular_velocity(centre, centre, orientation, half_length, cell, vec(0, 0.5, 0), R, repulsion, 2.0),
                   vec(0, 0, -8 * (0.1 - 2 * 0.5) / drag), TOL);
    // no stickiness: the motion of the cell does not matter, only the contact
    CHECK_VEC_NEAR(cell_angular_velocity(centre, centre, orientation, half_length, cell, vec(0, 0.5, 0), R, repulsion, 0.0),
                   vec(0, 0, -8 * 0.1 / drag), TOL);
    // no repulsion: only the traction
    CHECK_VEC_NEAR(cell_angular_velocity(centre, centre, orientation, half_length, cell, vec(0, 0.5, 0), R, 0.0, 2.0),
                   vec(0, 0, 8 * 2 * 0.5 / drag), TOL);
    // cell moving along the fibre: no traction torque
    CHECK_VEC_NEAR(cell_angular_velocity(centre, centre, orientation, half_length, cell, vec(0.5, 0, 0), R, repulsion, 2.0),
                   vec(0, 0, -8 * 0.1 / drag), TOL);
    // traction out of the plane (+z), contact in the plane
    CHECK_VEC_NEAR(cell_angular_velocity(centre, centre, orientation, half_length, cell, vec(0, 0, 0.5), R, repulsion, 2.0),
                   vec(0, -8 * 2 * 0.5 / drag, -8 * 0.1 / drag), TOL);
}

static void contact_torque_properties_in_3D()
{
    for (int n = 0; n < 1000; n++) {
        const Vec centre = random_vector(50), orientation = random_unit_vector();
        const double half_length = uniform(5, 40), R = uniform(6, 12), repulsion = uniform(0.1, 10);
        const Vec pivot = (n % 2 == 0) ? centre : add(centre, scale(uniform(-1, 1) * half_length, orientation));
        // a cell in contact with the fibre, anywhere along it (including beyond the ends)
        const double abscissa = uniform(-1.1, 1.1) * half_length;
        const Vec normal = random_perpendicular(orientation);
        const Vec cell = add(add(centre, scale(abscissa, orientation)), scale(uniform(0.05, 0.9) * R, normal));
        const Vec displacement = displacement_from_fibre(cell, centre, orientation, half_length);
        if (norm(displacement) >= R) {
            continue;
        }
        const Vec omega = cell_angular_velocity(pivot, centre, orientation, half_length, cell, vec(0, 0, 0), R, repulsion, 1.0);

        // perpendicular to the fibre
        CHECK_NEAR(dot(omega, orientation), 0.0, 1e-9);
        // the point of contact moves away from the cell
        const Vec impact = sub(cell, displacement);
        CHECK(dot(cross(omega, sub(impact, pivot)), displacement) <= 1e-12);

        // same result in another frame
        const Vec k = random_unit_vector(), shift = random_vector(100);
        const double angle = uniform(-PI, PI);
        const Vec moved = cell_angular_velocity(add(rotate(pivot, k, angle), shift), add(rotate(centre, k, angle), shift),
                                                rotate(orientation, k, angle), half_length, add(rotate(cell, k, angle), shift),
                                                vec(0, 0, 0), R, repulsion, 1.0);
        CHECK_VEC_NEAR(moved, rotate(omega, k, angle), 1e-9);
    }
}

// ---------------------------------------------------------------- main

struct Test
{
    const char* name;
    void (*function)();
};

#define TEST(function) { #function, function }

int main()
{
    const Test tests[] = {
        TEST(orientation_from_angles_closed_forms),
        TEST(orientation_from_angles_properties),
        TEST(nearest_point_along_the_fibre),
        TEST(nearest_point_beyond_the_ends),
        TEST(nearest_point_of_a_point_on_the_fibre),
        TEST(nearest_point_oblique_fibre),
        TEST(nearest_point_does_not_depend_on_the_sign_of_the_orientation),
        TEST(nearest_point_is_on_the_fibre_and_is_the_nearest),
        TEST(nearest_point_does_not_depend_on_the_frame),
        TEST(nearest_point_of_a_fibre_of_zero_length),
        TEST(angular_velocity_of_a_free_fibre),
        TEST(angular_velocity_of_a_crosslinked_fibre),
        TEST(angular_velocity_balances_the_drag_along_the_fibre),
        TEST(no_angular_velocity_without_torque),
        TEST(angular_velocities_of_several_cells_add_up),
        TEST(angular_velocity_properties_in_3D),
        TEST(rotation_closed_forms),
        TEST(rotation_keeps_the_sign_of_the_orientation),
        TEST(rotation_is_rigid),
        TEST(rotation_about_the_fibre_axis_changes_nothing),
        TEST(no_rotation_for_a_zero_rotation_vector),
        TEST(rotation_in_small_steps_equals_one_rotation),
        TEST(rotation_in_the_xy_plane_matches_the_2D_rotation),
        TEST(rotation_does_not_depend_on_the_frame),
        TEST(free_fibre_dragged_by_a_cell),
        TEST(crosslinked_fibre_dragged_by_a_cell),
        TEST(the_point_touched_by_the_cell_follows_the_cell),
        TEST(contact_velocity_closed_forms),
        TEST(contact_velocity_properties),
        TEST(cell_at_rest_pressing_on_a_fibre),
        TEST(cells_pressing_on_both_sides_of_a_fibre),
        TEST(contact_and_traction_add_up),
        TEST(contact_torque_properties_in_3D),
    };
    const int number_of_tests = sizeof(tests) / sizeof(tests[0]);

    int failed = 0;
    for (int i = 0; i < number_of_tests; i++) {
        checks_failed_in_test = 0;
        tests[i].function();
        if (checks_failed_in_test > 0) {
            failed++;
            std::printf("FAIL  %s (%d checks)\n", tests[i].name, checks_failed_in_test);
        }
        else {
            std::printf("ok    %s\n", tests[i].name);
        }
    }
    std::printf("\n%d tests, %d failed\n", number_of_tests, failed);
    return failed;
}
