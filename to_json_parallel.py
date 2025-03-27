import json
import multiprocessing as mp
import os
import time
import mmap
from collections import defaultdict

# Constants
NUM_PROCESSES = max(1, mp.cpu_count() - 1)  # Use all available cores except one

ATOM_TYPES = {}

with open("atom_info.json","r") as f:
    data = json.load(f)
    for i,each in enumerate(data):
        ATOM_TYPES[str(i+1)] = each

def process_chunk(chunk_data, atom_types_list, process_id):
    """Process a chunk of the file to extract atom data."""
    print(f"Process {process_id} starting on {len(chunk_data)} lines")
    start_time = time.time()
    
    result = defaultdict(list)
    current_snapshot = None
    current_timestep = None
    reading_atoms = False
    snapshot_count = 0
    
    lines = chunk_data.splitlines()
    line_count = len(lines)
    
    i = 0
    while i < line_count:
        line = lines[i].strip()
        
        if line.startswith("ITEM: TIMESTEP"):
            if current_snapshot:
                # End of snapshot, add to results
                for atom_type, atoms in current_snapshot.items():
                    # Sort atoms by ID
                    sorted_atoms = sorted(atoms, key=lambda x: int(x["atom_id"]))
                    result[atom_type].append(sorted_atoms)
                snapshot_count += 1
            
            reading_atoms = False
            current_snapshot = defaultdict(list)
            
            # Get timestep
            i += 1
            if i < line_count:
                current_timestep = int(lines[i].strip())
            
        elif line.startswith("ITEM: ATOMS"):
            reading_atoms = True
            
        elif reading_atoms and line and line[0].isdigit():
            parts = line.split()
            if len(parts) >= 5:
                atom_type = parts[1]
                if atom_type in atom_types_list:
                    atom_data = {
                        "atom_id": parts[0],
                        "atom_class": atom_type,
                        "atom_coordinate": [float(parts[2]), float(parts[3]), float(parts[4])]
                    }
                    current_snapshot[atom_type].append(atom_data)
        
        i += 1
    
    # Handle the last snapshot if any
    if current_snapshot:
        for atom_type, atoms in current_snapshot.items():
            sorted_atoms = sorted(atoms, key=lambda x: int(x["atom_id"]))
            result[atom_type].append(sorted_atoms)
        snapshot_count += 1
    
    end_time = time.time()
    print(f"Process {process_id} finished {snapshot_count} snapshots in {end_time - start_time:.2f} seconds")
    
    return result

def find_chunk_boundaries(file_path, num_chunks):
    """Find boundaries using memory mapping for faster processing."""
    file_size = os.path.getsize(file_path)
    
    with open(file_path, 'r+b') as f:
        # Memory-map the file for faster access
        mm = mmap.mmap(f.fileno(), 0)
        
        # Find all TIMESTEP positions
        timestep_positions = []
        pos = 0
        while True:
            pos = mm.find(b"ITEM: TIMESTEP", pos)
            if pos == -1:
                break
            timestep_positions.append(pos)
            pos += 1
        
        mm.close()
    
    # If we have fewer timesteps than chunks, adjust
    num_chunks = min(num_chunks, len(timestep_positions))
    if num_chunks <= 1:
        return [(0, file_size)]
    
    # Divide timesteps among chunks
    timesteps_per_chunk = len(timestep_positions) // num_chunks
    
    boundaries = []
    for i in range(num_chunks):
        start_idx = i * timesteps_per_chunk
        
        if i == num_chunks - 1:
            # Last chunk gets all remaining timesteps
            boundaries.append((timestep_positions[start_idx], file_size))
        else:
            end_idx = (i + 1) * timesteps_per_chunk
            boundaries.append((timestep_positions[start_idx], timestep_positions[end_idx]))
    
    return boundaries

def parallel_process_file(dump_file):
    """Process the entire file in parallel."""
    start_time = time.time()
    
    print(f"Starting parallel parsing with {NUM_PROCESSES} processes...")
    
    # Find chunk boundaries
    boundaries = find_chunk_boundaries(dump_file, NUM_PROCESSES)
    
    # Convert dict_keys to list for pickling
    atom_types_list = list(ATOM_TYPES.keys())
    
    # Read chunks in parallel
    tasks = []
    for i, (start, end) in enumerate(boundaries):
        print(f"Process {i+1} will handle bytes {start} to {end}")
        with open(dump_file, 'r') as f:
            f.seek(start)
            chunk_data = f.read(end - start)
        tasks.append((chunk_data, atom_types_list, i+1))
    
    # Set up process pool
    with mp.Pool(processes=NUM_PROCESSES) as pool:
        # Run processes
        results = pool.starmap(process_chunk, tasks)
    
    # Merge results
    merged_results = defaultdict(list)
    for result in results:
        for atom_type, snapshots in result.items():
            merged_results[atom_type].extend(snapshots)
    
    # Write results to JSON files
    for atom_type, atom_name in ATOM_TYPES.items():
        if atom_type in merged_results:
            output_file = f"{atom_name}.json"
            with open(output_file, 'w') as f:
                json.dump(merged_results[atom_type], f)
            print(f"Wrote {len(merged_results[atom_type])} snapshots to {output_file}")
    
    end_time = time.time()
    print(f"Parsing completed in {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    dump_file = 'out.dump'
    parallel_process_file(dump_file)
