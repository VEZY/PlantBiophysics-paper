"""
    RMSE(obs,sim)

Returns the Root Mean Squared Error between observations `obs` and simulations `sim`.
The closer to 0 the better.
"""
function RMSE(obs, sim, digits=2)
    return round(sqrt(sum((obs .- sim) .^ 2) / length(obs)), digits=digits)
end

"""
    nRMSE(obs,sim)

Returns the normalized Root Mean Squared Error between observations `obs` and simulations `sim`.
The closer to 0 the better.
"""
function nRMSE(obs, sim; digits=2)
    return round(
        sqrt(sum((obs .- sim) .^ 2) / length(obs)) / (findmax(obs)[1] - findmin(obs)[1]),
        digits=digits,
    )
end

"""
    EF(obs,sim)

Returns the Efficiency Factor between observations `obs` and simulations `sim` using NSE (Nash-Sutcliffe efficiency) model.
More information can be found at https://en.wikipedia.org/wiki/Nash%E2%80%93Sutcliffe_model_efficiency_coefficient.
The closer to 1 the better.
"""
function EF(obs, sim, digits=2)
    SSres = sum((obs - sim) .^ 2)
    SStot = sum((obs .- mean(obs)) .^ 2)
    return round(1 - SSres / SStot, digits=digits)
end

"""
	    Bias(obs,sim)

	Returns the bias between observations `obs` and simulations `sim`.
	The closer to 0 the better.
	"""
function Bias(obs, sim, digits=4)
    return round(mean(sim .- obs), digits=digits)
end

"""
	nBias(obs,sim; digits = 2)

Returns the normalised bias (%) between observations `obs` and simulations `sim`.
The closer to 0 the better.
"""
function nBias(obs, sim; digits=2)
    return round(mean((sim .- obs)) / (findmax(obs)[1] - findmin(obs)[1]), digits=digits)
end

"""
    load_medlyn_data(; type="curve", abs=0.85)

Load the Medlyn data using the working Figshare URL with fallback options.
Uses the PlantBiophysics.read_licor6400 function to parse the data.

# Arguments

- `type="curve"`: Specifies the type of data to load, either "curve" or "spot" measurements.
- `abs=0.85`: the leaf absorptance.
"""
function load_medlyn_data(; type="curve", abs=0.85)
    @assert type in ["curve", "spot"] "Type must be either 'curve' or 'spot'."
    file_id = type == "curve" ? "3402635" : "3402638"

    # Try the working Figshare URL first
    primary_url = "https://ndownloader.figshare.com/files/$file_id"

    # Keep alternative URLs as backup
    alternative_urls = [
        "https://figshare.com/ndownloader/files/$file_id",
        "https://figshare.com/download/file/$file_id"
    ]

    # Try primary URL first
    try
        println("Downloading from primary URL: $primary_url")
        return read_licor6400(Downloads.download(primary_url); abs=abs)
    catch e
        println("Failed with primary URL: $(e)")
    end

    # Try alternative URLs
    for url in alternative_urls
        try
            println("Trying alternative URL: $url")
            return read_licor6400(Downloads.download(url); abs=abs)
        catch e
            println("Failed with URL $url: $(e)")
            continue
        end
    end

    # If all downloads fail, provide helpful error message
    error("""
    Unable to download data from Figshare. This could be due to:
    1. Network connectivity issues
    2. Figshare server problems
    3. Changed file URLs

    Possible solutions:
    1. Check your internet connection
    2. Try again later
    3. Download the file manually from https://doi.org/10.6084/m9.figshare.1538079.v1
    """)
end
