import ee
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from PIL import Image
from tqdm import tqdm
import math

try:
    ee.Authenticate()
    ee.Initialize(project="earth-map-pixel")
except Exception as e:
    ee.Authenticate()
    ee.Initialize()
    #"Ukrain": [19.0, 40.0, 44.0, 40.0 ],

def generate_region_grid(region_name, preselected_regions, resolution=1000, custom_coords=None):

    if custom_coords is not None:
        region = ee.Geometry.Polygon(custom_coords)
        print(f"Custom region with coordinates: {custom_coords}")
    else:
        if region_name not in preselected_regions:
            raise ValueError(
                f"Region '{region_name}' not found. Available regions: {list(preselected_regions.keys())}")
        region = lon_min, lat_min, lon_max, lat_max = preselected_regions[region_name]
        region = ee.Geometry.Rectangle([lon_min, lat_min, lon_max, lat_max])
        print(f"Selected region '{region_name}' with coordinates: {lon_min, lat_min, lon_max, lat_max}")

    landcover = ee.Image('MODIS/006/MCD12Q1/2019_01_01').select('LC_Type1')
    landcover_region = landcover.clip(region)

    res_deg = resolution / 111320.0  # approximate conversion (111.32 km per degree)
    pixel_width = res_deg
    pixel_height = res_deg
    crs_transform = [pixel_width, 0, lon_min, 0, -pixel_height, lat_max]

    # Reproject the image using the new CRS transform.
    landcover_region = landcover_region.reproject(crs='EPSG:4326', crsTransform=crs_transform)
    sample_dict = landcover_region.sampleRectangle(region=region, defaultValue=0).getInfo()

    grid = np.array(sample_dict['properties']['LC_Type1'])
    return grid

def plot_grid(grid, color_map):

  # Create a list of colors in order of class 0 to 17.
  cmap_list = [color_map[i] for i in range(18)]
  cmap = mcolors.ListedColormap(cmap_list)

  plt.figure(figsize=(8, 6))
  # Use extent to map the array indices to geographic coordinates.
  plt.imshow(grid, cmap=cmap, interpolation='none')
  plt.title("MODIS Land Cover (2019)")
  plt.xlabel("Longitude")
  plt.ylabel("Latitude")
  plt.colorbar(ticks=range(18), label='Land Cover Class')

def create_img(grid, color_map):

  color_list = []
  for i in range(18):
      # Convert color name to an RGB tuple (values between 0 and 1)
      rgb = mcolors.to_rgb(color_map[i])
      # Scale to 0-255 and convert to integer.
      rgb = tuple(int(255 * x) for x in rgb)
      color_list.append(rgb)

  # ------------------------------
  # 7. Convert the grid to an RGB image.
  # Create an empty array for the RGB image.
  # ------------------------------
  height, width = grid.shape
  rgb_array = np.zeros((height, width, 3), dtype=np.uint8)

  # Map each pixel in the grid to its corresponding RGB color.
  for class_code, rgb in enumerate(color_list):
      mask = (grid == class_code)
      rgb_array[mask] = rgb

  # ------------------------------
  # 8. Save and display the image using Pillow.
  # ------------------------------
  img = Image.fromarray(rgb_array)
  img.save("lan)

def subdivide_region(lon_min, lat_min, lon_max, lat_max, n_tiles_x, n_tiles_y):
    """Divide the region into a list of subregions, each as [lon_min, lat_min, lon_max, lat_max]."""
    tile_width = (lon_max - lon_min) / n_tiles_x
    tile_height = (lat_max - lat_min) / n_tiles_y
    subregions = []
    for i in range(n_tiles_x):
        for j in range(n_tiles_y):
            sub_lon_min = lon_min + i * tile_width
            sub_lon_max = sub_lon_min + tile_width
            sub_lat_min = lat_min + j * tile_height
            sub_lat_max = sub_lat_min + tile_height
            subregions.append([sub_lon_min, sub_lat_min, sub_lon_max, sub_lat_max])
    return subregions

def global_coarse(region_name, resolution, preselected_regions, custom_region=None):

  if custom_region is not None:
    region = ee.Geometry.Polygon(custom_region)
    region_bounds = custom_region
    print(f"Custom region with coordinates: {custom_region}")
  else:
    if region_name not in preselected_regions:
      raise ValueError(f"Region '{region_name}' not found. Available regions: {list(preselected_regions.keys())}")
    region_bounds = lon_min, lat_min, lon_max, lat_max = preselected_regions[region_name]

  # Choose number of tiles (adjust as needed so that each tile’s pixel count is below the limit)
  n_tiles_x, n_tiles_y = 20, 10
  tiles = subdivide_region(*region_bounds, n_tiles_x, n_tiles_y)

  # Precompute your global image's CRS transform using your chosen resolution.
  res_deg = resolution / 111320.0
  pixel_width = res_deg
  pixel_height = res_deg
  # Use the global top-left corner (-180, 90) for the transform.
  crs_transform = [pixel_width, 0, -180, 0, -pixel_height, 90]

  # Reproject your global image once:
  landcover_global = landcover.clip(ee.Geometry.Rectangle(region_bounds)).reproject(
      crs='EPSG:4326', crsTransform=crs_transform
  )

  # Prepare a container for your global grid. Compute global grid dimensions:
  global_width = int(round((region_bounds[2] - region_bounds[0]) / pixel_width))
  global_height = int(round((region_bounds[3] - region_bounds[1]) / pixel_height))
  global_grid = np.empty((global_height, global_width), dtype=np.float32)

  # Process each tile:
  # with tqdm progress
  print("Number of tiles:", len(tiles))
  for i in tqdm(range(len(tiles))):

      tile = tiles[i]
      tile_geom = ee.Geometry.Rectangle(tile)
      sample_dict = landcover_global.sampleRectangle(
            region=tile_geom,
            defaultValue=0
        ).getInfo()
      # Assuming the band is nested under properties:
      tile_grid = np.array(sample_dict['properties']['LC_Type1'])

      # Compute expected dimensions for this tile:
      expected_tile_width = int(round((tile[2] - tile[0]) / pixel_width))
      expected_tile_height = int(round((tile[3] - tile[1]) / pixel_height))

      # Crop the tile grid to the expected dimensions if needed:
      tile_grid = tile_grid[:expected_tile_height, :expected_tile_width]

      # Determine the tile's position (offset) in the global grid.
      tile_lon_min, tile_lat_min, tile_lon_max, tile_lat_max = tile
      # Use math.floor for consistent rounding:
      col_start = int(math.floor((tile_lon_min - region_bounds[0]) / pixel_width))
      row_start = int(math.floor((region_bounds[3] - tile_lat_max) / pixel_height))

      # Determine insertion dimensions:
      h_tile, w_tile = tile_grid.shape
      available_h = global_grid.shape[0] - row_start
      available_w = global_grid.shape[1] - col_start
      insert_h = min(h_tile, available_h)
      insert_w = min(w_tile, available_w)

      if i == 0:
        print(np.unique(tile_grid))

      global_grid[row_start:row_start+insert_h, col_start:col_start+insert_w] = tile_grid[:insert_h, :insert_w]

  return global_grid

if __name__ == "__main__":
  preselected_regions = {
      "California": [-122.5, 37.0, -121.5, 38.0],
      "Berlin": [13.0, 52.3, 13.7, 52.7],
      "Sahara": [-10.0, 15.0, 30.0, 30.0],
      "earth": [-180.0, 180.0, -90.0, 90.0],
      "Ukrain": [19.0, 45.0, 44.0, 52.0 ],
      "NorthEast": [0.0, 0.0, 180.0, 90.0 ],
      "NorthWest": [-180.0, 0.0, 0.0, 90.0 ],
      "SouthWest": [-180.0, -90.0, 0.0, 0.0 ],
      "southEast": [0.0, -90.0, 180.0, 0.0 ],
      # Add more regions as needed...
  }
  color_mapping = {
      0: 'blue',             # Water
      1: 'darkgreen',        # Evergreen Needleleaf forest
      2: 'forestgreen',      # Evergreen Broadleaf forest
      3: 'limegreen',        # Deciduous Needleleaf forest
      4: 'lightgreen',       # Deciduous Broadleaf forest
      5: 'green',            # Mixed forest
      6: 'sienna',           # Closed shrublands
      7: 'peru',             # Open shrublands
      8: 'yellowgreen',      # Woody savannas
      9: 'gold',             # Savannas
      10: 'yellow',          # Grasslands
      11: 'cyan',            # Permanent wetlands
      12: 'orange',          # Croplands
      13: 'red',             # Urban and built-up
      14: 'pink',            # Cropland/Natural vegetation mosaic
      15: 'white',           # Snow and ice
      16: 'grey',            # Barren or sparsely vegetated
      17: 'black'            # Unclassified
  }

  selected_region = "southEast"
  resolution = 5000 # m per pixel, lower is more detailed, limit depends on region size
  #grid = generate_region_grid(selected_region, preselected_regions, resolution)
  grid = global_coarse(selected_region, resolution, preselected_regions)

  create_img(grid, color_mapping)
  plot_grid(grid, color_mapping) # depending on your pixel size this may fail

  np.savetxt('landcover_grid.csv', grid, delimiter=',', fmt='%d')
  print("Grid data saved to landcover_grid.csv")