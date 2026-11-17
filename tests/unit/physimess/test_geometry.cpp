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

// a fibre is symmetric: orientations u and -u describe the same fibre
static bool same_axis( const Vec& a, const Vec& b, double tolerance )
{
    return norm(cross(a, b)) <= tolerance && std::fabs(std::fabs(dot(a, b)) - 1.0) <= tolerance;
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

// angle of rotation applied by rotate_fibre_about_pivot, from its documentation
static double expected_angle( const Vec& lever, const Vec& force, double half_length, double stickiness )
{
    const double length = 2 * half_length;
    return stickiness * norm(cross(lever, force)) / (length * length / 12);
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

// ------------------------------------------------- rotation of a free fibre
// (the pivot is the centre of the fibre)

static void free_fibre_pushed_at_its_positive_end_in_plane()
{
    // fibre along x, force +y applied at x = +8: torque along +z,
    // anticlockwise rotation seen from +z
    const Vec centre = vec(0, 0, 0), impact = vec(8, 0, 0), force = vec(0, 0.5, 0);
    const double half_length = 10, stickiness = 1.0;
    const double angle = stickiness * 8 * 0.5 / (20.0 * 20.0 / 12); // 0.12 rad

    const FibreRotation rotation = rotate_fibre_about_pivot(centre, centre, impact, force, half_length, stickiness);
    CHECK(rotation.rotated);
    CHECK_VEC_NEAR(rotation.orientation, vec(std::cos(angle), std::sin(angle), 0), TOL);
    CHECK_VEC_NEAR(rotation.centre, centre, TOL);
}

static void free_fibre_pushed_at_its_negative_end_in_plane()
{
    // same force at x = -8: torque along -z, clockwise rotation seen from +z
    const Vec centre = vec(0, 0, 0), impact = vec(-8, 0, 0), force = vec(0, 0.5, 0);
    const double half_length = 10, stickiness = 1.0;
    const double angle = stickiness * 8 * 0.5 / (20.0 * 20.0 / 12);

    const FibreRotation rotation = rotate_fibre_about_pivot(centre, centre, impact, force, half_length, stickiness);
    CHECK(rotation.rotated);
    CHECK(same_axis(rotation.orientation, vec(std::cos(angle), -std::sin(angle), 0), TOL));
    // the returned orientation points from the pivot to the impact point
    CHECK_VEC_NEAR(rotation.orientation, vec(-std::cos(angle), std::sin(angle), 0), TOL);
    CHECK_VEC_NEAR(rotation.centre, centre, TOL);
}

static void free_fibre_pushed_out_of_plane()
{
    // fibre along x, force +z applied at x = +8: torque along -y,
    // the positive end of the fibre rises towards +z
    const Vec centre = vec(0, 0, 0), impact = vec(8, 0, 0), force = vec(0, 0, 0.5);
    const double half_length = 10, stickiness = 1.0;
    const double angle = stickiness * 8 * 0.5 / (20.0 * 20.0 / 12);

    FibreRotation rotation = rotate_fibre_about_pivot(centre, centre, impact, force, half_length, stickiness);
    CHECK(rotation.rotated);
    CHECK_VEC_NEAR(rotation.orientation, vec(std::cos(angle), 0, std::sin(angle)), TOL);

    // fibre along z, force +x applied at z = -8: the negative end moves towards +x
    rotation = rotate_fibre_about_pivot(centre, centre, vec(0, 0, -8), vec(0.5, 0, 0), half_length, stickiness);
    CHECK(rotation.rotated);
    CHECK(same_axis(rotation.orientation, vec(-std::sin(angle), 0, std::cos(angle)), TOL));
}

static void free_fibre_oblique_force()
{
    // the force (1, 2, 2) at x = +6 has the perpendicular part (0, 2, 2):
    // rotation of 6 * sqrt(8) * stickiness / (L^2/12) towards (0, 1, 1)/sqrt(2)
    const Vec centre = vec(0, 0, 0), impact = vec(6, 0, 0), force = vec(1, 2, 2);
    const double half_length = 12, stickiness = 0.8;
    const double angle = stickiness * 6 * std::sqrt(8.0) / (24.0 * 24.0 / 12);

    const FibreRotation rotation = rotate_fibre_about_pivot(centre, centre, impact, force, half_length, stickiness);
    CHECK(rotation.rotated);
    const double s = std::sin(angle) / std::sqrt(2.0);
    CHECK_VEC_NEAR(rotation.orientation, vec(std::cos(angle), s, s), TOL);
}

static void rotation_properties_in_3D()
{
    for (int n = 0; n < 1000; n++) {
        const Vec centre = random_vector(50), direction = random_unit_vector();
        const double half_length = uniform(5, 40), stickiness = uniform(0.1, 2);
        double abscissa = uniform(-half_length, half_length);
        if (std::fabs(abscissa) < 0.01) {
            abscissa = 0.01;
        }
        const Vec impact = add(centre, scale(abscissa, direction));
        const Vec lever = sub(impact, centre);
        // perpendicular and parallel parts of the force, chosen so that the angle is in (0.001, 0.5)
        const double angle = uniform(0.001, 0.5);
        const Vec perpendicular = random_perpendicular(direction);
        const double length = 2 * half_length;
        const double magnitude = angle * (length * length / 12) / (stickiness * std::fabs(abscissa));
        const Vec force = add(scale(magnitude, perpendicular), scale(uniform(-3, 3), direction));

        const FibreRotation rotation = rotate_fibre_about_pivot(centre, centre, impact, force, half_length, stickiness);
        CHECK(rotation.rotated);
        const Vec& orientation = rotation.orientation;

        // the orientation is a unit vector
        CHECK_NEAR(norm(orientation), 1.0, TOL);

        // the angle between the old and new axes is the documented one
        CHECK_NEAR(expected_angle(lever, force, half_length, stickiness), angle, 1e-10);
        const Vec old_direction = unit(lever);
        CHECK_NEAR(std::atan2(norm(cross(old_direction, orientation)), dot(old_direction, orientation)), angle, 1e-10);

        // the fibre rotates in the plane of the lever arm and the force (about the torque) ...
        const Vec torque = cross(lever, force);
        CHECK_NEAR(dot(orientation, unit(torque)), 0.0, 1e-10);
        // ... in the direction given by the right hand rule: the impact point follows the force
        CHECK(dot(cross(old_direction, orientation), torque) > 0);
        CHECK(dot(sub(orientation, old_direction), force) > 0);

        // it is the rotation by `angle` about the torque
        CHECK_VEC_NEAR(orientation, rotate(old_direction, unit(torque), angle), 1e-10);

        // a free fibre rotates about its centre, which does not move
        CHECK_VEC_NEAR(rotation.centre, centre, TOL);
    }
}

static void force_along_the_fibre_does_not_change_the_rotation()
{
    for (int n = 0; n < 200; n++) {
        const Vec centre = random_vector(50), direction = random_unit_vector();
        const double half_length = uniform(5, 40);
        const Vec impact = add(centre, scale(uniform(0.2, 1) * half_length, direction));
        const Vec force = scale(uniform(0.1, 1), random_perpendicular(direction));

        const FibreRotation reference = rotate_fibre_about_pivot(centre, centre, impact, force, half_length, 1.0);
        const FibreRotation other = rotate_fibre_about_pivot(centre, centre, impact,
                                                             add(force, scale(uniform(-5, 5), direction)), half_length, 1.0);
        CHECK(reference.rotated && other.rotated);
        CHECK_VEC_NEAR(other.orientation, reference.orientation, 1e-10);
    }
}

static void no_rotation_without_torque()
{
    const Vec centre = vec(1, 2, 3), direction = unit(vec(1, 2, -2));
    const double half_length = 10;
    const Vec impact = add(centre, scale(4, direction));

    // force applied at the pivot
    CHECK(!rotate_fibre_about_pivot(centre, centre, centre, vec(0, 1, 0), half_length, 1.0).rotated);
    // force along the fibre
    CHECK(!rotate_fibre_about_pivot(centre, centre, impact, scale(2, direction), half_length, 1.0).rotated);
    CHECK(!rotate_fibre_about_pivot(centre, centre, impact, scale(-2, direction), half_length, 1.0).rotated);
    // no force
    CHECK(!rotate_fibre_about_pivot(centre, centre, impact, vec(0, 0, 0), half_length, 1.0).rotated);
    // fibre which does not stick to the cell
    CHECK(!rotate_fibre_about_pivot(centre, centre, impact, vec(0, 1, 0), half_length, 0.0).rotated);
}

static void angle_is_proportional_to_stickiness_force_and_lever()
{
    const Vec centre = vec(0, 0, 0), force = vec(0, 0.3, 0.4); // norm 0.5, perpendicular to x
    const double half_length = 15; // L^2/12 = 75
    struct Case { double lever, stickiness, force_scale, angle; };
    const Case cases[] = {
        {5, 1.0, 1.0, 5 * 0.5 / 75.0},
        {10, 1.0, 1.0, 10 * 0.5 / 75.0},
        {5, 2.0, 1.0, 2 * 5 * 0.5 / 75.0},
        {5, 1.0, 3.0, 3 * 5 * 0.5 / 75.0},
        {-5, 0.5, 2.0, 5 * 0.5 / 75.0},
    };
    for (unsigned int i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        const Vec impact = vec(cases[i].lever, 0, 0);
        const FibreRotation rotation = rotate_fibre_about_pivot(centre, centre, impact,
                                                                scale(cases[i].force_scale, force), half_length, cases[i].stickiness);
        CHECK(rotation.rotated);
        const Vec old_direction = unit(impact);
        CHECK_NEAR(std::acos(dot(old_direction, rotation.orientation)), cases[i].angle, 1e-10);
    }
}

static void rotation_in_the_xy_plane_matches_the_2D_rotation()
{
    for (int n = 0; n < 500; n++) {
        const double phi = uniform(-PI, PI);
        const Vec centre = vec(uniform(-50, 50), uniform(-50, 50), 0);
        const Vec direction = vec(std::cos(phi), std::sin(phi), 0);
        const double half_length = uniform(5, 40), stickiness = uniform(0.1, 2);
        double abscissa = uniform(-half_length, half_length);
        if (std::fabs(abscissa) < 0.01) {
            abscissa = 0.01;
        }
        const Vec impact = add(centre, scale(abscissa, direction));
        const Vec force = vec(uniform(-0.5, 0.5), uniform(-0.5, 0.5), 0);
        const Vec lever = sub(impact, centre);

        const FibreRotation rotation = rotate_fibre_about_pivot(centre, centre, impact, force, half_length, stickiness);
        if (!rotation.rotated) {
            continue;
        }
        // 2D: signed angle from the z component of the torque, rotation matrix [c -s; s c]
        const double torque_z = lever[0] * force[1] - lever[1] * force[0];
        const double length = 2 * half_length;
        const double angle = stickiness * torque_z / (length * length / 12);
        const Vec old_direction = unit(lever);
        const Vec expected = vec(old_direction[0] * std::cos(angle) - old_direction[1] * std::sin(angle),
                                 old_direction[0] * std::sin(angle) + old_direction[1] * std::cos(angle), 0);

        CHECK_NEAR(rotation.orientation[2], 0.0, 0.0); // stays exactly in the plane
        CHECK_NEAR(rotation.centre[2], 0.0, 0.0);
        CHECK_VEC_NEAR(rotation.orientation, expected, 1e-10);
    }
}

static void rotation_does_not_depend_on_the_frame()
{
    for (int n = 0; n < 500; n++) {
        const Vec centre = random_vector(50), direction = random_unit_vector();
        const double half_length = uniform(5, 40), stickiness = uniform(0.1, 2);
        const Vec pivot = add(centre, scale(uniform(0, 1) * half_length, direction));
        const Vec impact = add(centre, scale(uniform(-1, 0) * half_length, direction));
        const Vec force = random_vector(0.3);
        const Vec k = random_unit_vector(), shift = random_vector(100);
        const double angle = uniform(-PI, PI);

        const FibreRotation reference = rotate_fibre_about_pivot(pivot, centre, impact, force, half_length, stickiness);
        const FibreRotation moved = rotate_fibre_about_pivot(add(rotate(pivot, k, angle), shift),
                                                             add(rotate(centre, k, angle), shift),
                                                             add(rotate(impact, k, angle), shift),
                                                             rotate(force, k, angle), half_length, stickiness);
        CHECK(reference.rotated == moved.rotated);
        if (reference.rotated && moved.rotated) {
            CHECK_VEC_NEAR(moved.orientation, rotate(reference.orientation, k, angle), 1e-9);
            CHECK_VEC_NEAR(moved.centre, add(rotate(reference.centre, k, angle), shift), 1e-9);
        }
    }
}

// ------------------------------------------- rotation of a crosslinked fibre
// (the pivot is the crosslink point, anywhere on the fibre)

// checks that the result is the rigid rotation of the fibre about the pivot
static void check_rigid_rotation_about_pivot( const Vec& pivot, const Vec& centre, const Vec& impact,
                                              const Vec& force, double half_length, double stickiness )
{
    const FibreRotation rotation = rotate_fibre_about_pivot(pivot, centre, impact, force, half_length, stickiness);
    CHECK(rotation.rotated);
    if (!rotation.rotated) {
        return;
    }
    const Vec lever = sub(impact, pivot);
    const Vec axis = unit(cross(lever, force));
    const double angle = expected_angle(lever, force, half_length, stickiness);

    // the fibre axis is rotated by `angle` about the torque
    CHECK_NEAR(norm(rotation.orientation), 1.0, TOL);
    CHECK_VEC_NEAR(rotation.orientation, rotate(unit(lever), axis, angle), 1e-10);

    // the centre is rotated about the pivot like the rest of the fibre
    CHECK_VEC_NEAR(rotation.centre, add(pivot, rotate(sub(centre, pivot), axis, angle)), 1e-9);
    // hence it stays at the same distance of the pivot, and the pivot stays on the fibre axis
    const Vec pivot_to_centre = sub(rotation.centre, pivot);
    CHECK_NEAR(norm(pivot_to_centre), norm(sub(centre, pivot)), 1e-9);
    CHECK_NEAR(norm(cross(pivot_to_centre, rotation.orientation)), 0.0, 1e-9);
}

static void fibre_hinged_at_one_end()
{
    // fibre from (0,0,0) to (20,0,0) crosslinked at the origin, pushed towards +y at x = 15
    const Vec pivot = vec(0, 0, 0), centre = vec(10, 0, 0), impact = vec(15, 0, 0), force = vec(0, 0.4, 0);
    const double half_length = 10, stickiness = 1.0;
    const double angle = 15 * 0.4 / (400.0 / 12); // 0.18 rad

    const FibreRotation rotation = rotate_fibre_about_pivot(pivot, centre, impact, force, half_length, stickiness);
    CHECK(rotation.rotated);
    CHECK_VEC_NEAR(rotation.orientation, vec(std::cos(angle), std::sin(angle), 0), TOL);
    CHECK_VEC_NEAR(rotation.centre, vec(10 * std::cos(angle), 10 * std::sin(angle), 0), TOL);

    // same fibre pushed towards +z: it leaves the plane, hinged at the origin
    const FibreRotation out_of_plane = rotate_fibre_about_pivot(pivot, centre, impact, vec(0, 0, 0.4), half_length, stickiness);
    CHECK(out_of_plane.rotated);
    CHECK_VEC_NEAR(out_of_plane.orientation, vec(std::cos(angle), 0, std::sin(angle)), TOL);
    CHECK_VEC_NEAR(out_of_plane.centre, vec(10 * std::cos(angle), 0, 10 * std::sin(angle)), TOL);
}

static void crosslinked_fibre_pushed_on_the_side_of_its_centre()
{
    for (int n = 0; n < 500; n++) {
        const Vec centre = random_vector(50), direction = random_unit_vector();
        const double half_length = uniform(5, 40);
        // pivot at abscissa p < 0, impact at abscissa q > p + 0.1 on the side of the centre
        const double p = uniform(-1, -0.05) * half_length;
        const double q = uniform(p + 0.1, half_length);
        const Vec pivot = add(centre, scale(p, direction)), impact = add(centre, scale(q, direction));
        const Vec force = add(scale(uniform(0.05, 0.3), random_perpendicular(direction)), scale(uniform(-1, 1), direction));
        check_rigid_rotation_about_pivot(pivot, centre, impact, force, half_length, uniform(0.1, 2));
    }
}

static void crosslinked_fibre_without_torque()
{
    const Vec pivot = vec(-4, 0, 0), centre = vec(0, 0, 0);
    // cell pushing exactly at the crosslink
    CHECK(!rotate_fibre_about_pivot(pivot, centre, pivot, vec(0, 1, 0), 10, 1.0).rotated);
    // cell pushing along the fibre
    CHECK(!rotate_fibre_about_pivot(pivot, centre, vec(5, 0, 0), vec(1, 0, 0), 10, 1.0).rotated);
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
        TEST(nearest_point_along_the_fibre),
        TEST(nearest_point_beyond_the_ends),
        TEST(nearest_point_of_a_point_on_the_fibre),
        TEST(nearest_point_oblique_fibre),
        TEST(nearest_point_does_not_depend_on_the_sign_of_the_orientation),
        TEST(nearest_point_is_on_the_fibre_and_is_the_nearest),
        TEST(nearest_point_does_not_depend_on_the_frame),
        TEST(free_fibre_pushed_at_its_positive_end_in_plane),
        TEST(free_fibre_pushed_at_its_negative_end_in_plane),
        TEST(free_fibre_pushed_out_of_plane),
        TEST(free_fibre_oblique_force),
        TEST(rotation_properties_in_3D),
        TEST(force_along_the_fibre_does_not_change_the_rotation),
        TEST(no_rotation_without_torque),
        TEST(angle_is_proportional_to_stickiness_force_and_lever),
        TEST(rotation_in_the_xy_plane_matches_the_2D_rotation),
        TEST(rotation_does_not_depend_on_the_frame),
        TEST(fibre_hinged_at_one_end),
        TEST(crosslinked_fibre_pushed_on_the_side_of_its_centre),
        TEST(crosslinked_fibre_without_torque),
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
