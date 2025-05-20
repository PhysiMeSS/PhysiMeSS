# Merged PhysiCell to ParaView Converter
# Based on work by Othmane Hayoun-Mya (2025-04-04)
#
# This script converts PhysiCell simulation output files (.mat and .xml)
# into ParaView-compatible formats (.vtu and .pvd) for visualization.

import os
import glob
import scipy.io
import numpy as np
import xml.etree.ElementTree as ET
from vtk import (vtkUnstructuredGrid, vtkPoints, vtkVertex, vtkDoubleArray,
                vtkCellArray, vtkXMLUnstructuredGridWriter, vtkPolyData,
                vtkSphereSource, vtkCylinderSource, vtkAppendPolyData,
                vtkTransform, vtkTransformPolyDataFilter)
import re

def parse_physiboss_labels(xml_file):
    """Parse the XML file to get labels and metadata for each field"""
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        
        # Find the labels section in the XML
        labels_elem = root.find(".//labels")
        if labels_elem is None:
            raise ValueError(f"Could not find labels in XML file: {xml_file}")
        
        # Parse each label
        labels = {}
        for label in labels_elem.findall('label'):
            index = int(label.get('index'))
            size = int(label.get('size'))
            units = label.get('units', 'none')
            name = label.text.strip() if label.text else f"field_{index}"
            
            # Sanitize name to ensure XML compatibility
            name = ''.join(c if c.isalnum() or c in '_- ' else '_' for c in name)
            
            labels[index] = {'name': name, 'size': size, 'units': units}
        
        return labels
    
    except Exception as e:
        print(f"Error parsing XML file {xml_file}: {str(e)}")
        return {}


def mat_to_vtu(mat_file, xml_file, output_file):
    """Convert a PhysiCell .mat file to VTU format with proper labels"""
    try:
        # Get labels from XML
        labels = parse_physiboss_labels(xml_file)
        if not labels:
            return False
        
        # Load mat file
        data = scipy.io.loadmat(mat_file)
        cells_data = data['cells']
        
        # Debug print
        print(f"Processing {os.path.basename(mat_file)}")
        print(f"Cells array shape: {cells_data.shape}")
        
        # Get number of cells
        num_cells = cells_data.shape[1]
        
        # Create VTK grid
        grid = vtkUnstructuredGrid()
        points = vtkPoints()
        
        # Track valid cells
        valid_cell_indices = []
        cell_types = []  # Store cell types for later
        cell_radii = []  # Store cell radii for later
        
        # Process each cell
        for i in range(num_cells):
            try:
                # Get cell position
                x = float(cells_data[1][i])
                y = float(cells_data[2][i])
                z = float(cells_data[3][i])
                
                # Get cell properties (if available)
                cell_type = int(cells_data[5][i]) if cells_data.shape[0] > 5 else 0
                radius = float(cells_data[4][i]) if cells_data.shape[0] > 4 else 10.0
                
                # Add point
                points.InsertNextPoint(x, y, z)
                
                # Track this cell
                valid_cell_indices.append(i)
                cell_types.append(cell_type)
                cell_radii.append(radius)
                
            except (IndexError, ValueError) as e:
                print(f"Warning: Error processing cell {i}: {e}")
                continue
        
        grid.SetPoints(points)
        
        # Create vertices
        vertices = vtkCellArray()
        for i in range(len(valid_cell_indices)):
            vertex = vtkVertex()
            vertex.GetPointIds().SetId(0, i)
            vertices.InsertNextCell(vertex)
        grid.SetCells(vtkVertex().GetCellType(), vertices)
        
        # Add cell type and radius as point data
        type_array = vtkDoubleArray()
        type_array.SetName("cell_type")
        type_array.SetNumberOfComponents(1)
        
        radius_array = vtkDoubleArray()
        radius_array.SetName("radius")
        radius_array.SetNumberOfComponents(1)
        
        for i in range(len(valid_cell_indices)):
            type_array.InsertNextValue(cell_types[i])
            radius_array.InsertNextValue(cell_radii[i])
        
        grid.GetPointData().AddArray(type_array)
        grid.GetPointData().AddArray(radius_array)
        
        # Add other cell data
        for idx, label_info in labels.items():
            name = label_info['name']
            size = label_info['size']
            
            # Skip position data as it's already handled
            if name == "position":
                continue
                
            # Skip indices that don't exist in the original data
            if idx >= cells_data.shape[0]:
                continue
                
            try:
                # Create array for this field
                array = vtkDoubleArray()
                array.SetName(name)
                
                if size == 1:
                    # Scalar data
                    array.SetNumberOfComponents(1)
                    for cell_idx in valid_cell_indices:
                        try:
                            value = float(cells_data[idx][cell_idx])
                            if np.isnan(value) or np.isinf(value):
                                value = 0.0
                            array.InsertNextValue(value)
                        except (IndexError, ValueError):
                            array.InsertNextValue(0.0)
                else:
                    # Vector data
                    array.SetNumberOfComponents(size)
                    for cell_idx in valid_cell_indices:
                        try:
                            vector = []
                            for j in range(size):
                                if idx + j < cells_data.shape[0]:
                                    val = float(cells_data[idx + j][cell_idx])
                                    if np.isnan(val) or np.isinf(val):
                                        val = 0.0
                                    vector.append(val)
                                else:
                                    vector.append(0.0)
                            array.InsertNextTuple(vector)
                        except (IndexError, ValueError):
                            array.InsertNextTuple([0.0] * size)
                
                grid.GetPointData().AddArray(array)
            except Exception as e:
                print(f"Warning: Skipping field {name} due to error: {str(e)}")
                continue
        
        # Write to file
        writer = vtkXMLUnstructuredGridWriter()
        writer.SetFileName(output_file)
        writer.SetInputData(grid)
        
        # Use binary format and compression for better file handling
        writer.SetDataModeToBinary()
        writer.SetCompressorTypeToZLib()
        
        writer.Write()
        print(f"Successfully wrote {os.path.basename(output_file)}")
        return True
        
    except Exception as e:
        print(f"Error processing file {os.path.basename(mat_file)}: {str(e)}")
        return False

def create_pvd(output_dir, clean=False, file_prefix="timestep"):
    """Create a ParaView Data (PVD) file indexing all timesteps"""
    # Clean existing output if requested
    if clean:
        # Remove existing VTU files with the specified prefix
        vtu_files = glob.glob(os.path.join(output_dir, f"{file_prefix}_*.vtu"))
        for file in vtu_files:
            print(f"Removing existing VTU file: {file}")
            os.remove(file)
        
        # Remove existing PVD file
        pvd_file = os.path.join(output_dir, "simulation.pvd")
        if os.path.exists(pvd_file):
            print(f"Removing existing PVD file: {pvd_file}")
            os.remove(pvd_file)
    
    # Find all mat files
    mat_files = sorted(glob.glob(os.path.join(output_dir, "*_cells.mat")))
    if not mat_files:
        raise ValueError(f"No .mat files found in {output_dir}")
    
    # PVD file header
    pvd_content = ['<?xml version="1.0"?>',
                   '<VTKFile type="Collection" version="0.1">',
                   '  <Collection>']
    
    successful_conversions = 0
    vtu_files = []
    
    # Process each timestep
    for i, mat_file in enumerate(mat_files):
        # Get corresponding XML file
        xml_file = mat_file.replace("_cells.mat", ".xml")
        if not os.path.exists(xml_file):
            print(f"Warning: XML file not found for {mat_file}")
            continue
        
        # Extract timestep from filename
        base_name = os.path.basename(mat_file)
        match = re.search(r'(\d+)_cells\.mat$', base_name)
        if not match:
            print(f"Warning: Could not extract timestep from filename: {base_name}")
            continue
            
        timestep = int(match.group(1))
        
        # Create VTU filename directly in the output directory
        vtu_file = os.path.join(output_dir, f"{file_prefix}_{timestep:06d}.vtu")
        
        print(f"\nProcessing timestep {timestep}...")
        
        # Convert mat to vtu
        if mat_to_vtu(mat_file, xml_file, vtu_file):
            # Use just the basename in the PVD file
            vtu_basename = os.path.basename(vtu_file)
            pvd_content.append(f'    <DataSet timestep="{timestep}" group="" part="0" file="{vtu_basename}"/>')
            successful_conversions += 1
            vtu_files.append(vtu_file)
    
    # Close PVD file
    pvd_content.extend(['  </Collection>',
                       '</VTKFile>'])
    
    # Write PVD file only if we have successful conversions
    if successful_conversions > 0:
        pvd_file = os.path.join(output_dir, "simulation.pvd")
        with open(pvd_file, 'w') as f:
            f.write('\n'.join(pvd_content))
        print(f"\nSuccessfully converted {successful_conversions} timesteps")
        print(f"Output files saved to:")
        print(f"  - PVD file: {os.path.abspath(pvd_file)}")
        print(f"  - VTU files: {len(vtu_files)} files in {os.path.abspath(output_dir)}")
        return pvd_file, vtu_files
    else:
        print("\nNo successful conversions")
        return None, []

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Convert PhysiCell output to ParaView format')
    parser.add_argument('output_dir', help='Directory containing .mat and .xml output files')
    parser.add_argument('--clean', action='store_true', help='Remove existing output files before processing')
    parser.add_argument('--prefix', default='timestep', help='File prefix for VTU files (default: timestep)')
    
    args = parser.parse_args()
    
    pvd_file, vtu_files = create_pvd(args.output_dir, args.clean, args.prefix)
    
    if pvd_file:
        print(f"\nConversion complete! To visualize these files in ParaView:")
        print(f"1. Open ParaView")
        print(f"2. File > Open > Navigate to: {os.path.abspath(pvd_file)}")
        print(f"3. Click 'Apply' in the Properties panel to load the data")

if __name__ == '__main__':
    main()