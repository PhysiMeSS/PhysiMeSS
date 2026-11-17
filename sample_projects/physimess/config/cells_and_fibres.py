# ParaView script: cells as spheres and fibres as cylinders.
#
# The simulation output must first be converted with paraview_viewer.py
# ("make vtu output=<output folder>"), which writes simulation.pvd and one
# .vtu file per saved time step in the output folder.
#
# In ParaView (graphical interface):
#   1. File > Open, select simulation.pvd in the output folder, click Apply
#   2. View > Python Shell, click "Run Script" and select this file
#
# From the command line, with the pvpython program installed with ParaView:
#   pvpython cells_and_fibres.py <output folder>/simulation.pvd
#       --png image.png       save an image of the last time step (or of --step N)
#       --frames <folder>     save one image per time step (frame_0000.png, ...)
#
# Cells are spheres of their own radius. Fibres are cylinders along their
# orientation, with the length and radius of their cell definition
# (fibre_length, fibre_radius). Both are coloured by cell type.
#
# Tested with ParaView 6.2.

import os
import sys

from paraview.simple import *

FIBRE_SIDES = 16         # number of sides of the cylinders
SPHERE_RESOLUTION = 24   # resolution of the spheres


def find_input(arguments):
    """The data to display: the file given on the command line, or the source selected in ParaView."""
    files = [a for a in arguments if a.lower().endswith((".pvd", ".vtu"))]
    if files:
        if files[0].lower().endswith(".pvd"):
            return PVDReader(FileName=os.path.abspath(files[0]))
        return XMLUnstructuredGridReader(FileName=[os.path.abspath(f) for f in files])
    source = GetActiveSource()
    if source is None:
        raise RuntimeError("Open simulation.pvd first (File > Open, then Apply), or give it on the command line.")
    return source


def option(arguments, name):
    """Value following the option `name` on the command line, or None."""
    if name in arguments and arguments.index(name) + 1 < len(arguments):
        return arguments[arguments.index(name) + 1]
    return None


def array_range(source, name):
    """(min, max) of a point data array over the current time step, or None if it does not exist."""
    info = source.GetPointDataInformation().GetArray(name)
    return info.GetRange(0) if info is not None else None


def select(source, array, value, name):
    """The points for which `array` equals `value`."""
    threshold = Threshold(Input=source, registrationName=name)
    threshold.Scalars = ["POINTS", array]
    threshold.LowerThreshold = value - 0.5
    threshold.UpperThreshold = value + 0.5
    return threshold


def colour_by_type(display, view):
    ColorBy(display, ("POINTS", "cell_type"))
    display.SetScalarBarVisibility(view, False)


def show_cells_and_fibres(source, view):
    source.UpdatePipeline()
    if array_range(source, "cell_type") is None or array_range(source, "radius") is None:
        raise RuntimeError("The selected source has no cell data: select simulation.pvd in the pipeline browser.")
    Hide(source, view)

    has_fibre_flag = array_range(source, "is_fibre") is not None
    if not has_fibre_flag:
        print("No 'is_fibre' data (files converted with an older paraview_viewer.py): all agents are drawn as spheres.")

    # cells: spheres
    cells = select(source, "is_fibre", 0, "cells") if has_fibre_flag else source
    spheres = Glyph(Input=cells, GlyphType="Sphere", registrationName="spheres")
    spheres.GlyphType.Radius = 1.0
    spheres.GlyphType.ThetaResolution = SPHERE_RESOLUTION
    spheres.GlyphType.PhiResolution = SPHERE_RESOLUTION
    spheres.OrientationArray = ["POINTS", "No orientation array"]
    spheres.ScaleArray = ["POINTS", "radius"]
    spheres.ScaleFactor = 1.0
    spheres.GlyphMode = "All Points"
    colour_by_type(Show(spheres, view), view)

    # fibres: segments along the orientation, of length fibre_length, thickened into cylinders
    if has_fibre_flag and array_range(source, "is_fibre")[1] > 0.5:
        fibres = select(source, "is_fibre", 1, "fibres")
        segments = Glyph(Input=fibres, GlyphType="Line", registrationName="fibre_axes")
        segments.OrientationArray = ["POINTS", "orientation"]
        segments.ScaleArray = ["POINTS", "fibre_length"]
        segments.ScaleFactor = 1.0
        segments.GlyphMode = "All Points"
        cylinders = Tube(Input=segments, registrationName="cylinders")
        radius = array_range(source, "fibre_radius")
        cylinders.Radius = radius[1] if radius is not None and radius[1] > 0 else 2.0
        cylinders.NumberofSides = FIBRE_SIDES
        cylinders.Capping = 1
        colour_by_type(Show(cylinders, view), view)

    view.ResetCamera()
    Render(view)


def main(arguments):
    source = find_input(arguments)
    view = GetActiveViewOrCreate("RenderView")
    show_cells_and_fibres(source, view)

    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    times = list(scene.TimeKeeper.TimestepValues) if scene.TimeKeeper.TimestepValues else [0.0]

    png = option(arguments, "--png")
    frames = option(arguments, "--frames")
    if png or frames:
        view.ViewSize = [1200, 900]
    if png:
        step = option(arguments, "--step")
        scene.AnimationTime = times[int(step) if step is not None else -1]
        view.ResetCamera()
        SaveScreenshot(os.path.abspath(png), view)
        print("Saved", os.path.abspath(png))
    if frames:
        os.makedirs(frames, exist_ok=True)
        scene.AnimationTime = times[-1]
        view.ResetCamera()   # same camera for all the frames, framed on the last time step
        for i, time in enumerate(times):
            scene.AnimationTime = time
            SaveScreenshot(os.path.join(os.path.abspath(frames), "frame_%04d.png" % i), view)
        print("Saved", len(times), "images in", os.path.abspath(frames))


main(sys.argv[1:])
