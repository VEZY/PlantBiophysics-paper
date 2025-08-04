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
    load_medlyn_data()

Load the Medlyn data using the working Figshare URL with fallback options.
Uses the PlantBiophysics.read_licor6400 function to parse the data.
"""
function load_medlyn_data()
    # Try the working Figshare URL first
    primary_url = "https://ndownloader.figshare.com/files/3402635"

    # Keep alternative URLs as backup
    alternative_urls = [
        "https://figshare.com/ndownloader/files/3402635",
        "https://figshare.com/download/file/3402635"
    ]

    # Try primary URL first
    try
        println("Downloading from primary URL: $primary_url")
        return read_licor6400(Downloads.download(primary_url))
    catch e
        println("Failed with primary URL: $(e)")
    end

    # Try alternative URLs
    for url in alternative_urls
        try
            println("Trying alternative URL: $url")
            return read_licor6400(Downloads.download(url))
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
    3. Download the file manually from https://figshare.com/articles/dataset/Tumbarumba_Gas_Exchange/1538079
    """)
end
