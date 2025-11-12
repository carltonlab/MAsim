import numpy as np
from scipy import interpolate

def load_map_data(filename):
    """Load the mapping data from a file."""
    return np.loadtxt(filename)

class GenomicPositionGenerator:
    def __init__(self, map_data):
        """Initialize the generator with map data."""
        self.genomic_positions, self.genetic_positions = map_data[:, 0], map_data[:, 1]
        
        # Create the interpolation function once during initialization
        self.interp_func = interpolate.interp1d(self.genetic_positions, self.genomic_positions, 
                                                kind='linear', fill_value='extrapolate')
        
        # Store the maximum genetic position for efficient random number generation
        self.max_genetic_pos = self.genetic_positions[-1]

    def generate_position(self):
        """Generate a random genomic position based on genetic distance."""
        # Generate a random genetic position
        random_genetic_pos = np.random.random() * self.max_genetic_pos
        
        # Use the interpolation function to get the corresponding genomic position
        return float(self.interp_func(random_genetic_pos))

# Usage
#map_data = load_map_data('ivs.txt')
#generator = GenomicPositionGenerator(map_data)

# Generate a single random position
#random_genomic_position = generator.generate_position()
#print(f"Random genomic position: {random_genomic_position}")

# Generate multiple positions (more efficiently now)
#num_positions = 1000
#random_positions = [generator.generate_position() for _ in range(num_positions)]

# If you need to generate a large number of positions at once, you can further optimize:
#many_random_genetic_pos = np.random.random(num_positions) * generator.max_genetic_pos
#many_random_genomic_pos = generator.interp_func(many_random_genetic_pos)
