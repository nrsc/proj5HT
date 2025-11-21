import numpy as np
import matplotlib.pyplot as plt
import pyreadr
import pynwb
from pathlib import Path
import rdata
import os
import h5py
import argparse
import sys


def read_nwb_traces(nwb_fname, trace_nums, chan_num=0):
    """
    Read NWB traces for specified trace numbers.
    
    Parameters:
    -----------
    nwb_fname : str or Path
        Path to the NWB file
    trace_nums : list
        List of trace numbers to read
    chan_num : int
        Channel number (default: 0)
        
    Returns:
    --------
    dict : Dictionary containing stim, acq, time, and ttl data
    """
    
    # Open the MIES file
    h5_file = h5py.File(nwb_fname, 'r')
    # Location of the stimulus and acquisition data within the NWB/MIES file
    stim_path = 'stimulus/presentation'
    stim_keys = list(h5_file[stim_path].keys())
    try:
        acq_path = 'acquisition/timeseries'
        acq_keys = list(h5_file[acq_path].keys())
    except:
        acq_path = 'acquisition'
        acq_keys = list(h5_file[acq_path].keys())
    
    stim = []
    acq  = []
    ttl  = []
    time = []
    
    # Cycle through the selected traces and add them to the list
    for trace_ct, trace_num in enumerate(trace_nums):
        # Grab the stimulus and recording traces
        this_stim_key = 'data_%s_DA%i'     % ('{0:05d}'.format(trace_num), chan_num)
        this_ttl_key  = 'data_%s_TTL1_1'  % ('{0:05d}'.format(trace_num))
        this_acq_key  = 'data_%s_AD%i'  % ('{0:05d}'.format(trace_num), chan_num)
        stim_trace_num = stim_keys.index(this_stim_key) if this_stim_key in stim_keys else -1
        ttl_trace_num  = stim_keys.index(this_ttl_key)  if this_ttl_key  in stim_keys else -1
        acq_trace_num  = acq_keys.index(this_acq_key)   if this_acq_key  in acq_keys  else -1
        stim.append(np.array(h5_file[stim_path + '/' + stim_keys[stim_trace_num] + '/data'])) if stim_trace_num >= 0 else stim.append(np.empty(0))
        ttl.append(np.array(h5_file[stim_path   + '/' + stim_keys[ttl_trace_num]  + '/data'])) if ttl_trace_num  >= 0 else ttl.append(np.empty(0))
        acq.append(np.array(h5_file[acq_path   + '/' + acq_keys[acq_trace_num]   + '/data'])) if acq_trace_num  >= 0 else acq.append(np.empty(0))
        
        # Identify the sampling interval
        trace_si = 1 / h5_file[acq_path + '/' + acq_keys[acq_trace_num] + '/starting_time'].attrs['rate']
        
        # Generate an array of time stamps
        time.append(np.arange(0, len(stim[-1])) * trace_si)
        
    # Close the file
    h5_file.close()
    
    return {'stim': stim, 'acq': acq, 'time': time, 'ttl': ttl}


def generate_cell_figure(cell_name, root_path=None, output_path=None):
    """
    Generate and save a figure for a specific cell.
    
    Parameters:
    -----------
    cell_name : str
        Name of the cell to process
    root_path : str or Path, optional
        Root directory path (default: //allen/programs/celltypes/workgroups/hct/SawchukS/rookery)
    output_path : str or Path, optional
        Directory to save the figure (default: root_path / cell_name)
        
    Returns:
    --------
    str : Path to the saved figure
    """
    
    # Set default root path if not provided
    if root_path is None:
        root = Path("//allen/programs/celltypes/workgroups/hct/SawchukS/rookery")
    else:
        root = Path(root_path)
    
    # Set default output path if not provided
    if output_path is None:
        output_dir = root / cell_name
    else:
        output_dir = Path(output_path)
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Processing cell: {cell_name}")
    print(f"Root directory: {root}")
    print(f"Output directory: {output_dir}")
    
    # Construct file path
    pattern = f"{cell_name}-srt.rds"
    file_path = root / cell_name / pattern
    
    # Check if RDS file exists
    if not file_path.exists():
        raise FileNotFoundError(f"RDS file does not exist: {file_path}")
    
    # Load RDS data
    print("Loading RDS data...")
    srt = rdata.read_rds(str(file_path))
    print("Loaded RDS data successfully!")
    
    # Get NWB path
    nwb_path = Path(srt["files"]["nwb"][0])
    
    # Get spike sweeps
    spike_sweeps = srt["dfs"]["spikeTTL"]["spike_TTL_sweep_numbers"].astype(int) - 1
    
    # Read NWB traces
    print("Reading NWB traces...")
    tst = read_nwb_traces(nwb_path, spike_sweeps)
    tst0 = np.concatenate(tst["acq"])
    
    # Prepare time data
    time0 = []
    offset = srt["dfs"]["spikeTTL"]["spike_puff_output"]["time"][1]
    
    for t in tst["time"]:
        time0.append(t + offset)
        offset += t[-1]  # add last value of current vector as offset
    
    time0 = np.concatenate(time0)
    
    # Get additional data for subplots
    x2 = srt["dfs"]["spikeTTL"]["spike_puff_output"]["time"]
    y2 = srt["dfs"]["spikeTTL"]["spike_puff_output"]["instRate"]
    y3 = srt["dfs"]["spikeTTL"]["spike_puff_output"]["percent_change"]
    
    # Create the figure
    print("Generating figure...")
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
    
    # Top plot - voltage trace
    ax1.plot(time0, tst0)
    ax1.set_ylabel('mV')
    ax1.set_title(f'{cell_name} - Voltage Trace and Analysis')
    ax1.grid(True, alpha=0.3)
    
    # Middle plot - instantaneous frequency
    ax2.plot(x2, y2, color='red')
    ax2.set_ylabel('Instantaneous Frequency (Hz)')
    ax2.grid(True, alpha=0.3)
    
    # Bottom plot - percent change
    ax3.plot(x2, y3, color='blue')
    ax3.set_xlabel('Time (seconds)')
    ax3.set_ylabel('% Change')
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save the figure
    output_file = output_dir / f"{cell_name}_analysis.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Figure saved to: {output_file}")
    
    # Also save as PDF for better quality
    output_file_pdf = output_dir / f"{cell_name}_analysis.pdf"
    plt.savefig(output_file_pdf, bbox_inches='tight')
    print(f"PDF version saved to: {output_file_pdf}")
    
    plt.close()  # Close the figure to free memory
    
    return str(output_file)


def main():
    """Main function to handle command line arguments."""
    parser = argparse.ArgumentParser(description='Generate analysis figure for a specific cell.')
    parser.add_argument('cell_name', help='Name of the cell to process')
    parser.add_argument('--root', '-r', help='Root directory path', default=None)
    parser.add_argument('--output', '-o', help='Output directory path', default=None)
    parser.add_argument('--show', '-s', action='store_true', help='Display the figure after saving')
    
    args = parser.parse_args()
    
    try:
        output_file = generate_cell_figure(
            cell_name=args.cell_name,
            root_path=args.root,
            output_path=args.output
        )
        
        if args.show:
            # Load and display the saved figure
            from PIL import Image
            img = Image.open(output_file)
            img.show()
            
        print(f"Successfully processed {args.cell_name}")
        
    except Exception as e:
        print(f"Error processing {args.cell_name}: {e}")
        sys.exit(1)


if __name__ == "__main__":
    # If running as script with command line arguments
    if len(sys.argv) > 1:
        main()
    else:
        # Example usage if no arguments provided
        print("Usage examples:")
        print("python build_loadData.py QN24.26.001.4A.17.02")
        print("python build_loadData.py QN24.26.001.4A.17.02 --root /custom/root/path")
        print("python build_loadData.py QN24.26.001.4A.17.02 --output /custom/output/path")
        print("python build_loadData.py QN24.26.001.4A.17.02 --show")
        print("\nOr use as a module:")
        print("from build_loadData import generate_cell_figure")
        print("generate_cell_figure('QN24.26.001.4A.17.02')")





