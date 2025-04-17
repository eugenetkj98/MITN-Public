"""
Author: Eugene Tan
Date Created: 14/4/2024
Last Updated: 14/4/2024
Function to lookup and linearly interpolate value from a raster at a given latitude and longitude
"""

module RasterLookup

export interp_lookup

using GeoIO
using Rasters

"""
    interp_lookup(raster::Raster, lat::Float64, lon::Float64)

Lookup and linearly (simplex) interpolate for a value of a raster at the given latitude and longitude
"""

function interp_lookup(raster, lat, lon)

    # Get raster reference positions
    raster_lats = lookup(raster, Y)
    raster_lons = lookup(raster, X)

    # Lookup raster idx reference location
    lat_idx = argmin(abs.(raster_lats .- lat))
    lon_idx = argmin(abs.(raster_lons .- lon))

    # Get reference raster value
    ref_val = raster[lon_idx,lat_idx]
    ref_lat = raster_lats[lat_idx]
    ref_lon = raster_lons[lon_idx]

    # Get fiducial values for interpolation on triangle
    δlat = lat-ref_lat
    δlon = lon-ref_lon

    neighbour_lat = NaN
    neighbour_lat_val = NaN
    if δlat > 0 # ref is south of obs
        neighbour_lat = raster_lats[lat_idx-1]
        neighbour_lat_val = raster[lon_idx, lat_idx-1]
    else
        neighbour_lat = raster_lats[lat_idx+1]
        neighbour_lat_val = raster[lon_idx, lat_idx+1]
    end

    neighbour_lon = NaN
    neighbour_lon_val = NaN
    if δlon <0
        neighbour_lon = raster_lons[lon_idx-1]
        neighbour_lon_val = raster[lon_idx-1,lat_idx]
    else
        neighbour_lon = raster_lons[lon_idx+1]
        neighbour_lon_val = raster[lon_idx+1,lat_idx]
    end

    # Perform interpolation
    Δlat = neighbour_lat - ref_lat
    J = δlat/Δlat
    Δlon = neighbour_lon - ref_lon
    I = δlon/Δlon

    interp_val = (1-I-J)*ref_val + I*neighbour_lon_val + J*neighbour_lat_val

    return interp_val
end

end