
#SA environment
import os
import geopandas as gpd
import matplotlib.pyplot as plt
import networkx as nx
import osmnx as ox
from shapely.geometry import LineString
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.geometry import Point, LineString, Polygon


from scipy.spatial import ConvexHull
from shapely.geometry import box
ox.__version__


import pyrosm
from pyrosm import OSM, get_data

import osmnet

import pandana as pdna

from pandana.loaders import osm

import geopandas as gpd
import pandas as pd


city="Manama_Prj"

prem_outputs=r"C:/Users/Kochulem/Dropbox (Personal)/Projects/Data/"

#Kabul
#Bus_path = r'D:\OneDrive - United Nations\03 Geospatial Data\Edwin\POS_Transport_Degurba_extent_data_ready_to_process\Afghanistan\Kabul_Updated_Data\Kabul_Bus_Stops.shp'
#Road_path = r'D:\OneDrive - United Nations\03 Geospatial Data\Edwin\POS_Transport_Degurba_extent_data_ready_to_process\Afghanistan\Kabul_Updated_Data\Kabul_Roads.shp'


#Manama
Bus_path = r'D:\OneDrive - United Nations\03 Geospatial Data\Edwin\POS_Transport_Degurba_extent_data_ready_to_process\Bahrain\Manama\Manama_combined_Updated_Bustops.shp'  # Replace 'path/to/shapefile.shp' with the actual path
Road_path = r'D:\OneDrive - United Nations\03 Geospatial Data\Edwin\POS_Transport_Degurba_extent_data_ready_to_process\Bahrain\Manama\Manama_Roads.shp'  # Replace 'path/to/shapefile.shp' with the actual path


#Tirane
#Road_path = r'D:\OneDrive - United Nations\03 Geospatial Data\Edwin\Processed_Statistics\Processed_stats_11.2.1_and_11.7.1_final\Albania\Tirana\Roads_1.shp'
#Bus_path = r'D:\OneDrive - United Nations\03 Geospatial Data\Edwin\Processed_Statistics\Processed_stats_11.2.1_and_11.7.1_final\Albania\Tirana\Bustops_1.shp'


#Gaborone
#Bus_path = r'D:\OneDrive - United Nations\03 Geospatial Data\Edwin\Processed_Statistics\Processed_stats_11.2.1_and_11.7.1_final\Botswana\Gaborone\Bustops_1.shp'
#Road_path = r'D:\OneDrive - United Nations\03 Geospatial Data\Edwin\Processed_Statistics\Processed_stats_11.2.1_and_11.7.1_final\Botswana\Gaborone\Roads_1.shp'


output_geopackage = r"C:/Users/Kochulem/Dropbox (Personal)/Projects/Data/"+city+"_isochrone_three.gpkg"

output_SHP=r"C:/Users/Kochulem/Dropbox (Personal)/Projects/Data/"+city+"_isochrone_Bus_SA.shp"

output_Point_SHP = r"C:/Users/Kochulem/Dropbox (Personal)/Projects/Data/"+city+"_Points.shp"

output_Line_SHP=r"C:/Users/Kochulem/Dropbox (Personal)/Projects/Data/"+city+"_Lines.shp"

#Boundary of region of interest.
Road_gdf=gpd.read_file(Road_path)

current_crs=Road_gdf.crs


# Extract the EPSG code from the current CRS
epsg_code = int(current_crs.to_epsg())

print("EPSG code:", epsg_code)

Road_gdf_geog = Road_gdf.to_crs(epsg=4326)





left, bottom, right, top = Road_gdf_geog.total_bounds
bbox = top, bottom, right, left
poly = ox.utils_geo.bbox_to_poly(*bbox)

G = ox.graph_from_polygon(poly, network_type='all')


# Project the graph network to a particular EPSG
#G = ox.project_graph(G, to_crs='EPSG:'+str(epsg_code))




print(G)

# Define travel parameters
trip_distances = [500]  # in meters (500 meters-access to bus stops, 1000meters-access to high capacity stops)
travel_speed = 4.8  # walking speed in km/hour

snap_tolerance = 100  # Adjust this value as needed


# Calculate meters per minute for travel speed
meters_per_minute = travel_speed * 1000 / 60  # km per hour to meters per minute

# Add an edge attribute for distance in meters required to traverse each edge
for _, _, _, data in G.edges(data=True, keys=True):
    data["distance"] = data["length"]  # Assuming the length attribute represents the edge length in meters

# Read the shapefile into a GeoDataFrame
gdf_Bus2 = gpd.read_file(Bus_path)

# Project the GeoDataFrame to EPSG 4326 (if needed)
gdf_Bus = gdf_Bus2.to_crs(epsg=4326)





# Identify nearest network nodes for each location of interest
for index, location in gdf_Bus.iterrows():
    # Snap facility locations to the nearest nodes in the graph network
    nearest_node = ox.distance.nearest_nodes(G, location.geometry.x, location.geometry.y)
    gdf_Bus.loc[index, "nearest_node"] = nearest_node




# By now, you have already identified nearest nodes and stored them in gdf_Bus
# Initialize an empty list to store the central nodes for isochrones
central_nodes = gdf_Bus["nearest_node"].tolist()

# Initialize an empty list to store isochrone polygons
isochrone_polys = []

# Loop through each central node and generate isochrones
for central_node in central_nodes:
    # Initialize an empty list to store node points for each isochrone
    node_points = []
    # Generate subgraph for the current central node and specified travel distances
    for trip_distance in sorted(trip_distances, reverse=True):
        subgraph = nx.ego_graph(G, central_node, radius=trip_distance, distance="distance")

        # Extract node coordinates and create Point geometries
        node_points.extend([Point((data["x"], data["y"])) for node, data in subgraph.nodes(data=True)])

        # Create convex hull from node points and add to isochrone_polys list
        bounding_poly = gpd.GeoSeries(node_points).unary_union.convex_hull
        isochrone_polys.append(bounding_poly)



print(isochrone_polys)


# Convert the Polygon to a GeoDataFrame
gdf_isochrones = gpd.GeoDataFrame(geometry=isochrone_polys)


# Filter out only the polygons
polygons_gdf = gdf_isochrones[gdf_isochrones["geometry"].geom_type == "Polygon"]


Point_gdf = gdf_isochrones[gdf_isochrones["geometry"].geom_type == "Point"]

print(polygons_gdf)


Lines_gdf = gdf_isochrones[gdf_isochrones["geometry"].geom_type == "LineString"]

# Check if there are any polygons in the GeoDataFrame
if len(polygons_gdf) > 0:
    # Define the coordinate reference system (CRS) for the GeoDataFrame
    ##polygons_gdf.crs = 'EPSG:'+str(epsg_code)  # Set CRS to WGS84 (latitude/longitude)
    polygons_gdf.crs = 'EPSG:4326'


    Projected_polygons_gdf = polygons_gdf.to_crs(epsg=epsg_code)

    # Export the GeoDataFrame containing polygons to a shapefile
    output_shapefile = output_SHP
    Projected_polygons_gdf.to_file(output_shapefile)
    print("Shapefile containing only polygons exported successfully.")
else:
    print("No polygons found in the GeoDataFrame.")




# Check if there are any polygons in the GeoDataFrame
if len(Point_gdf) > 0:
    # Define the coordinate reference system (CRS) for the GeoDataFrame
    ##Point_gdf.crs = 'EPSG:'+str(epsg_code)  # Set CRS to WGS84 (latitude/longitude)
    Point_gdf.crs = 'EPSG:4326'


    Projected_points_gdf = Point_gdf.to_crs(epsg=epsg_code)

    # Export the GeoDataFrame containing polygons to a shapefile
    output_shapefile_Point = output_Point_SHP
    Projected_points_gdf.to_file(output_shapefile_Point)
    print("Shapefile containing only points exported successfully.")
else:
    print("No points found in the GeoDataFrame.")




# Check if there are any polygons in the GeoDataFrame
if len(Lines_gdf) > 0:
    # Define the coordinate reference system (CRS) for the GeoDataFrame
    ##Lines_gdf.crs = 'EPSG:'+str(epsg_code)  # Set CRS to WGS84 (latitude/longitude)
    Lines_gdf.crs = 'EPSG:4326'


    Projected_Lines_gdf = Lines_gdf.to_crs(epsg=epsg_code)

    # Export the GeoDataFrame containing polygons to a shapefile
    output_shapefile_Line = output_Line_SHP
    Projected_Lines_gdf.to_file(output_shapefile_Line)
    print("Shapefile containing only Lines exported successfully.")
else:
    print("No Lines found in the GeoDataFrame.")




#End of the old code that was working