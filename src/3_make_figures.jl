# Generate maps of species richness, predicted numbers of links (median and 68% PI) and stability score

## Map of species richness

# Function to read species richness from a dataset in .tif extension
# Usable with data from biodiversitymapping.org
function read_richness(data)
    # Read dataset in .tif extension using ArchGDAL
    raw_tif = ArchGDAL.read(ArchGDAL.read(data))

    # Convert in matrix format and transpose
    raw = Matrix(transpose(Float32.(raw_tif[:, :, 1])))
    raw = raw[end:-1:1, :]

    # We are working on the actual pixels, which results in a big map, but we will
    # coarsen it later
    sp = SimpleSDMPredictor(raw, -16671202.0, 16818798.0, -8418528.0, 8381472.0)

    # We want to remove the 255 (0xff) and 65535 (oxffff)(NA values)
    sp_255 = filter(x -> sp.grid[x] in [255, 65535], CartesianIndices(sp.grid))

    for i in LinearIndices(sp.grid)[sp_255]
        sp.grid[i] = 0
    end

    # Return raster
    return(sp)
end

# Data paths
# Species richness of terrestrial mammals, birds and amphibians
# Rasters obtained from QGIS after uniformizing their extents
mammals_path = "data/BiodiversityMapping/mammals.tif"
birds_path = "data/BiodiversityMapping/birds.tif"
amphibians_path = "data/BiodiversityMapping/amphibians.tif"

# Read species richness of specific taxons (terrestrial mammals, birds and amphibians)
sp_mammals = read_richness(mammals_path)
sp_birds = read_richness(birds_path)
sp_amphibians = read_richness(amphibians_path)

# Raster sizes
size(sp_mammals)
size(sp_birds)
size(sp_amphibians)

# Visualize species richness
heatmap(sp_mammals)
heatmap(sp_birds)
heatmap(sp_amphibians)

# New SimpleSDMResponse with same coordinates
total_sp_richness = similar(sp_birds)

# Sum of birds, mammals and amphibians
for i in LinearIndices(total_sp_richness.grid)
    total_sp_richness.grid[i] = sp_mammals.grid[i] + sp_birds.grid[i] + sp_amphibians.grid[i]
end

# Change 0 to NAs
sp_NA = filter(x -> total_sp_richness.grid[x] == 0, CartesianIndices(total_sp_richness.grid))

for i in LinearIndices(total_sp_richness.grid)[sp_NA]
    total_sp_richness.grid[i] = NaN
end

# Map
# TK: bounding box at the beginning of this script
xlim = (-10500000, -3800000)
ylim = (5000000, 8381472)

heatmap(total_sp_richness, dpi=300, clim=(0,400), xlim=xlim, ylim=ylim, frame=:box)
savefig("figures/species_richness.png")


## Maps of predicted numbers of links (median and 68% PI) and stability score

# Load measures
measures = CSV.read(joinpath("data", "measures.csv"), DataFrame)

# Where are there at least 5 species?
cidx = filter(x -> total_sp_richness.grid[x] >= 5, CartesianIndices(total_sp_richness.grid))
cidx0 = filter(x -> total_sp_richness.grid[x] < 5, CartesianIndices(total_sp_richness.grid))

# SimpleSDMResponse objects
L_median_SDM = similar(total_sp_richness)
L_PI_SDM = similar(total_sp_richness)
stab_SDM = similar(total_sp_richness)

# Change the values in the SimpleSDMResponse objects by their respective measures according to species richness
function map_measures(measure, SDMObject)
    for i in LinearIndices(total_sp_richness.grid)[cidx]
    SDMObject.grid[i] = measures[measures.S .== total_sp_richness.grid[i], measure][1]
    end
    for i in LinearIndices(total_sp_richness.grid)[cidx0]
        SDMObject.grid[i] = NaN
    end
end

measures[!, :L_median_log] .= log10.(measures.L_median)
measures[!, :L_PI_log] .= log10.(measures.L_PI)

map_measures("L_median_log", L_median_SDM)
map_measures("L_PI_log", L_PI_SDM)
map_measures("stab", stab_SDM)

# Maps

heatmap(L_median_SDM, color=:BrBG_11, clim=(0,4), xlim=xlim, ylim=ylim, dpi=300, framestyle=:box) # log scale
savefig(joinpath("figures", "maps_L_median.png"))

heatmap(L_PI_SDM, color=:BrBG_11, clim=(0,4), xlim=xlim, ylim=ylim, dpi=300, framestyle=:box) # log scale
savefig(joinpath("figures", "maps_L_PI.png"))

heatmap(stab_SDM, color=:RdBu_11, dpi=300,  xlim=xlim, ylim=ylim, framestyle=:box)
savefig(joinpath("figures", "maps_stab.png"))
