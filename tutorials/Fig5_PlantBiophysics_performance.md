~~~
<!-- PlutoStaticHTML.Begin -->
<!--
    # This information is used for caching.
    [PlutoStaticHTML.State]
    input_sha = "eae3e5b41aedd65b40a0022cca12967b61b8888bbb3148ac2f483339bdac9fa1"
    julia_version = "1.8.5"
-->

<div class="markdown"><h1>PlantBiophysics.jl benchmark</h1>
<p>The main objective of this notebook is to compare the computational times of <code>PlantBiophysics.jl</code> against the <a href="https://github.com/RemkoDuursma/plantecophys">plantecophys</a> R package and the <a href="https://github.com/cropbox/LeafGasExchange.jl">LeafGasExchange.jl</a> Julia package from the <a href="https://github.com/cropbox/Cropbox.jl">Cropbox.jl</a> framework. The comparison follows three steps:</p>
<ul>
<li><p>create an N-large basis of random conditions.</p>
</li>
<li><p>benchmark the computational time of the three packages via similar functions &#40;<em>i.e.</em> photosynthesis-stomatal conductance-energy balance coupled model for C3 leaves&#41;: <code>energy_balance</code>, <code>photosynEB</code> and <code>simulate</code> with <code>ModelC3MD</code>.</p>
</li>
<li><p>compare the results with plots and statistics.</p>
</li>
</ul>
<p>This notebook does not perform the benchmark by itself for the obvious reason that it takes forever to run, and because there is an overhead cause by Pluto &#40;benchmark and reactive are not a good mix&#41;. Instead, it shows the outputs of <a href="https://github.com/VEZY/PlantBiophysics-paper/blob/main/tutorials/Fig5_PlantBiophysics_performance_noPluto.jl">a script from this repository</a> that implements the code shown here. If you want to perform the benchmark by yourself, you can run this script from the command line. The versions used for the above dependencies are available in the <a href="https://github.com/VEZY/PlantBiophysics-paper/blob/main/tutorials/Project.toml">Project.toml</a> of the repository.</p>
<h2>Importing the dependencies:</h2>
<h6>Note</h6>
<p>Make sure to have R installed on your computer first.</p>
<p>Loading the Julia packages:</p>
<pre><code class="language-julia">begin
    using CairoMakie
	using BenchmarkTools
	using PlantBiophysics
    using Cropbox
    using LeafGasExchange
    using RCall
end</code></pre>
</div>

<pre class='language-julia'><code class='language-julia'>begin
    using Statistics
    using DataFrames
    using CSV
    using Random
    using PlantBiophysics
end</code></pre>



<div class="markdown"><h2>Parameters</h2>
<h3>Benchmark parameters</h3>
<p>You&#39;ll find below the main parameters of the benchmark. In a few words, each package runs a simulation for <code>N</code> different time-steps <code>microbenchmark_steps</code> times repeated <code>microbenchmark_evals</code> times. We make <code>N</code> different simulations because the simulation duration can vary depending on the inputs due to iterative computations in the code, <em>i.e.</em> different initial conditions can make the algorithms converge more or less rapidly.</p>
</div>

<pre class='language-julia'><code class='language-julia'>begin
    Random.seed!(1) # Set random seed
    microbenchmark_steps = 100 # Number of times the microbenchmark is run
    microbenchmark_evals = 1 # N. times each sample is run to be sure of the output
    N = 100 # Number of timesteps simulated for each microbenchmark step
end</code></pre>
<pre id='var-N' class='code-output documenter-example-output'>100</pre>


<div class="markdown"><h3>Random input simulation dataset</h3>
<p>We create possible ranges for input parameters. These ranges where chosen so all of the three packages don&#39;t return errors during computation &#40;<code>plantecophys</code> has issues with low temperatures&#41;.</p>
<ul>
<li><p>Ta: air temperature &#40;<span class="tex">$°C$</span>&#41;</p>
</li>
<li><p>Wind: wind speed &#40;<span class="tex">$m.s^&#123;-1&#125;$</span>&#41;</p>
</li>
<li><p>P: ambient pressure &#40;<span class="tex">$kPa$</span>&#41;</p>
</li>
<li><p>Rh: relative humidity &#40;between 0 and 1&#41;</p>
</li>
<li><p>Ca: air CO₂ concentration &#40;<span class="tex">$ppm$</span>&#41;</p>
</li>
<li><p>Jmax: potential rate of electron transport &#40;<span class="tex">$\mu mol_&#123;CO2&#125;.m^&#123;-2&#125;.s^&#123;-1&#125;$</span>&#41;</p>
</li>
<li><p>Vmax: maximum rate of Rubisco activity &#40;<span class="tex">$\mu mol_&#123;CO2&#125;.m^&#123;-2&#125;.s^&#123;-1&#125;$</span>&#41;</p>
</li>
<li><p>Rd: mitochondrial respiration in the light at reference temperature &#40;<span class="tex">$\mu mol_&#123;CO2&#125;.m^&#123;-2&#125;.s^&#123;-1&#125;$</span>&#41;</p>
</li>
<li><p>TPU: triose phosphate utilization-limited photosynthesis rate &#40;<span class="tex">$\mu mol_&#123;CO2&#125;.m^&#123;-2&#125;.s^&#123;-1&#125;$</span>&#41;</p>
</li>
<li><p>Rs: short-wave net radiation &#40;<span class="tex">$W.m^&#123;-1&#125;$</span>&#41;</p>
</li>
<li><p>skyF: Sun-visible fraction of the leaf &#40;between 0 and 1&#41;</p>
</li>
<li><p>d: characteristic length &#40;<span class="tex">$m$</span>&#41;</p>
</li>
<li><p>g0: residual stomatal conductance &#40;<span class="tex">$mol_&#123;CO2&#125;.m^&#123;-2&#125;.s^&#123;-1&#125;$</span>&#41;</p>
</li>
<li><p>g1: slope of the stomatal conductance relationship.</p>
</li>
</ul>
</div>

<pre class='language-julia'><code class='language-julia'>begin
    # Create the ranges of input parameters
    length_range = 10000
    Rs = range(10, 500, length=length_range)
    Ta = range(18, 40, length=length_range)
    Wind = range(0.5, 20, length=length_range)
    P = range(90, 101, length=length_range)
    Rh = range(0.1, 0.98, length=length_range)
    Ca = range(360, 900, length=length_range)
    skyF = range(0.0, 1.0, length=length_range)
    d = range(0.001, 0.5, length=length_range)
    Jmax = range(200.0, 300.0, length=length_range)
    Vmax = range(150.0, 250.0, length=length_range)
    Rd = range(0.3, 2.0, length=length_range)
    TPU = range(5.0, 20.0, length=length_range)
    g0 = range(0.001, 2.0, length=length_range)
    g1 = range(0.5, 15.0, length=length_range)
    vars = hcat([Ta, Wind, P, Rh, Ca, Jmax, Vmax, Rd, Rs, skyF, d, TPU, g0, g1])
    nothing
end</code></pre>



<div class="markdown"><p>We then sample <code>N</code> conditions from the given ranges:</p>
</div>

<pre class='language-julia'><code class='language-julia'>begin
    set = [rand.(vars) for i = 1:N]
    set = reshape(vcat(set...), (length(set[1]), length(set)))'
    name = [
        "T",
        "Wind",
        "P",
        "Rh",
        "Ca",
        "JMaxRef",
        "VcMaxRef",
        "RdRef",
        "Rs",
        "sky_fraction",
        "d",
        "TPURef",
        "g0",
        "g1",
    ]
    set = DataFrame(set, name)
    @. set[!, :vpd] = e_sat(set.T) - vapor_pressure(set.T, set.Rh)
    @. set[!, :PPFD] = set.Rs * 0.48 * 4.57
    set
end</code></pre>
<table>
<tr>
<th></th>
<th>T</th>
<th>Wind</th>
<th>P</th>
<th>Rh</th>
<th>Ca</th>
<th>JMaxRef</th>
<th>VcMaxRef</th>
<th>RdRef</th>
<th>...</th>
</tr>
<tr>
<td>1</td>
<td>39.8174</td>
<td>2.17132</td>
<td>98.4026</td>
<td>0.599626</td>
<td>644.824</td>
<td>226.203</td>
<td>246.79</td>
<td>1.20194</td>
<td></td>
</tr>
<tr>
<td>2</td>
<td>27.6458</td>
<td>5.4847</td>
<td>94.8724</td>
<td>0.143212</td>
<td>765.689</td>
<td>282.758</td>
<td>226.498</td>
<td>0.458796</td>
<td></td>
</tr>
<tr>
<td>3</td>
<td>26.363</td>
<td>19.7289</td>
<td>96.6898</td>
<td>0.487151</td>
<td>506.193</td>
<td>214.641</td>
<td>176.273</td>
<td>0.671147</td>
<td></td>
</tr>
<tr>
<td>4</td>
<td>29.7228</td>
<td>2.37414</td>
<td>96.7602</td>
<td>0.500176</td>
<td>850.153</td>
<td>205.851</td>
<td>196.945</td>
<td>0.347945</td>
<td></td>
</tr>
<tr>
<td>5</td>
<td>30.9703</td>
<td>11.698</td>
<td>93.2365</td>
<td>0.935908</td>
<td>581.854</td>
<td>286.579</td>
<td>202.385</td>
<td>1.98674</td>
<td></td>
</tr>
<tr>
<td>6</td>
<td>23.791</td>
<td>1.988</td>
<td>98.0726</td>
<td>0.628317</td>
<td>520.828</td>
<td>222.662</td>
<td>182.513</td>
<td>0.372087</td>
<td></td>
</tr>
<tr>
<td>7</td>
<td>35.0341</td>
<td>11.6512</td>
<td>99.6282</td>
<td>0.190561</td>
<td>462.286</td>
<td>219.232</td>
<td>179.783</td>
<td>0.755306</td>
<td></td>
</tr>
<tr>
<td>8</td>
<td>36.0308</td>
<td>10.6196</td>
<td>91.9142</td>
<td>0.354081</td>
<td>841.998</td>
<td>291.259</td>
<td>245.8</td>
<td>1.95393</td>
<td></td>
</tr>
<tr>
<td>9</td>
<td>38.0858</td>
<td>5.67387</td>
<td>100.636</td>
<td>0.109417</td>
<td>395.428</td>
<td>256.396</td>
<td>175.663</td>
<td>0.788799</td>
<td></td>
</tr>
<tr>
<td>10</td>
<td>23.5402</td>
<td>4.19952</td>
<td>96.5578</td>
<td>0.146381</td>
<td>429.559</td>
<td>244.084</td>
<td>175.703</td>
<td>1.39559</td>
<td></td>
</tr>
<tr>
<td>...</td>
</tr>
<tr>
<td>100</td>
<td>23.3443</td>
<td>6.82838</td>
<td>97.7492</td>
<td>0.612211</td>
<td>545.563</td>
<td>276.058</td>
<td>198.445</td>
<td>1.86518</td>
<td></td>
</tr>
</table>



<div class="markdown"><h2>Benchmarking</h2>
</div>


<div class="markdown"><h5>plantecophys</h5>
<p>Preparing R to make the benchmark:</p>
<pre><code class="language-julia">R&quot;&quot;&quot;
if&#40;&#33;require&#40;&quot;plantecophys&quot;&#41;&#41;&#123;
    install.packages&#40;&quot;plantecophys&quot;, repos &#61; &quot;https://cloud.r-project.org&quot;&#41;
&#125;
if&#40;&#33;require&#40;&quot;microbenchmark&quot;&#41;&#41;&#123;
    install.packages&#40;&quot;microbenchmark&quot;, repos &#61; &quot;https://cloud.r-project.org&quot;&#41;
&#125;
&quot;&quot;&quot;

# Make variables available to the R session
@rput set N microbenchmark_steps
</code></pre>
<p>Making the benchmark:</p>
<pre><code class="language-julia">R&quot;&quot;&quot;
# Define the function call in a function that takes a list as input to limit DataFrame overhead
function_EB &lt;- function&#40;input&#41; &#123;
    PhotosynEB&#40;
        Tair &#61; input&#36;Tair, VPD &#61; input&#36;VPD, Wind &#61; input&#36;Wind,
        Wleaf &#61; input&#36;Wleaf,Ca &#61; input&#36;Ca,  StomatalRatio &#61; 1,
        LeafAbs &#61; input&#36;LeafAbs, gsmodel &#61; &quot;BBOpti&quot;, g0 &#61; input&#36;g0, g1 &#61; input&#36;g1,
        alpha &#61; 0.24, theta &#61; 0.7, Jmax &#61; input&#36;Jmax,
        Vcmax &#61; input&#36;Vcmax, TPU &#61; input&#36;TPU, Rd &#61; input&#36;Rd,
        RH &#61; input&#36;RH, PPFD&#61;input&#36;PPFD, Patm &#61; input&#36;Patm
    &#41;
&#125;

time_PE &#61; c&#40;&#41;
for&#40;i in seq_len&#40;N&#41;&#41;&#123;
    # Put the inputs into a vector to limit dataframe overhead:
    input &#61; list&#40;
        Tair &#61; set&#36;T&#91;i&#93;, VPD &#61; set&#36;vpd&#91;i&#93;, Wind &#61; set&#36;Wind&#91;i&#93;, Wleaf &#61; set&#36;d&#91;i&#93;,
        Ca &#61; set&#36;Ca&#91;i&#93;, LeafAbs &#61; set&#36;sky_fraction&#91;i&#93;, g0 &#61; set&#36;g0&#91;i&#93;, g1 &#61; set&#36;g1&#91;i&#93;,
        Jmax &#61; set&#36;JMaxRef&#91;i&#93;, Vcmax &#61; set&#36;VcMaxRef&#91;i&#93;, TPU &#61; set&#36;TPURef&#91;i&#93;,
        Rd &#61; set&#36;RdRef&#91;i&#93;, RH &#61; set&#36;Rh&#91;i&#93;*100, PPFD&#61;set&#36;PPFD&#91;i&#93;,Patm &#61; set&#36;P&#91;i&#93;
    &#41;

    m &#61; microbenchmark&#40;function_EB&#40;input&#41;, times &#61; microbenchmark_steps&#41;

    time_PE &#61; append&#40;time_PE,m&#36;time * 10e-9&#41; # transform in seconds
&#125;
&quot;&quot;&quot;

@rget time_PE</code></pre>
</div>


<div class="markdown"><h5>LeafGasExchange.jl</h5>
<p>Note that we benchmark <code>LeafGasExchange.jl</code> with the <code>nounit</code> flag to compute a fair comparison with <code>PlantBiophysics.jl</code> in case computing units takes time &#40;it shouldn&#39;t much&#41;.</p>
<pre><code class="language-julia">time_LG &#61; &#91;&#93;
n_lg &#61; fill&#40;0, N&#41;
for i &#61; 1:N
    config &#61;
        :Weather &#61;&gt; &#40;
            PFD&#61;set.PPFD&#91;i&#93;,
            CO2&#61;set.Ca&#91;i&#93;,
            RH&#61;set.Rh&#91;i&#93; * 100,
            T_air&#61;set.T&#91;i&#93;,
            wind&#61;set.Wind&#91;i&#93;,
            P_air&#61;set.P&#91;i&#93;,
            g0&#61;set.g0&#91;i&#93;,
            g1&#61;set.g1&#91;i&#93;,
            Vcmax&#61;set.VcMaxRef&#91;i&#93;,
            Jmax&#61;set.JMaxRef&#91;i&#93;,
            Rd&#61;set.RdRef&#91;i&#93;,
            TPU&#61;set.TPURef&#91;i&#93;,
        &#41;
    b_LG &#61;
        @benchmark simulate&#40;&#36;ModelC3MD; config&#61;&#36;config&#41; evals &#61; microbenchmark_evals samples &#61;
            microbenchmark_steps
    append&#33;&#40;time_LG, b_LG.times .* 1e-9&#41; # transform in seconds
    n_lg&#91;i&#93; &#61; 1
end</code></pre>
</div>


<div class="markdown"><h5>PlantBiophysics.jl</h5>
<p>Benchmarking <code>PlantBiophysics.jl</code>:</p>
<pre><code class="language-julia">constants &#61; Constants&#40;&#41;
time_PB &#61; &#91;&#93;
for i &#61; 1:N
    leaf &#61; ModelList&#40;
        energy_balance&#61;Monteith&#40;&#41;,
        photosynthesis&#61;Fvcb&#40;
            VcMaxRef&#61;set.VcMaxRef&#91;i&#93;,
            JMaxRef&#61;set.JMaxRef&#91;i&#93;,
            RdRef&#61;set.RdRef&#91;i&#93;,
            TPURef&#61;set.TPURef&#91;i&#93;,
        &#41;,
        stomatal_conductance&#61;Medlyn&#40;set.g0&#91;i&#93;, set.g1&#91;i&#93;&#41;,
        status&#61;&#40;
            Rₛ&#61;set.Rs&#91;i&#93;,
            sky_fraction&#61;set.sky_fraction&#91;i&#93;,
            PPFD&#61;set.PPFD&#91;i&#93;,
            d&#61;set.d&#91;i&#93;,
        &#41;,
    &#41;
    deps &#61; PlantSimEngine.dep&#40;leaf&#41;
    meteo &#61; Atmosphere&#40;T&#61;set.T&#91;i&#93;, Wind&#61;set.Wind&#91;i&#93;, P&#61;set.P&#91;i&#93;, Rh&#61;set.Rh&#91;i&#93;, Cₐ&#61;set.Ca&#91;i&#93;&#41;
    st &#61; PlantMeteo.row_struct&#40;leaf.status&#91;1&#93;&#41;
    b_PB &#61; @benchmark run&#33;&#40;&#36;leaf, &#36;deps, &#36;st, &#36;meteo, &#36;constants, nothing&#41; evals &#61;
        microbenchmark_evals samples &#61; microbenchmark_steps
    append&#33;&#40;time_PB, b_PB.times .* 1e-9&#41; # transform in seconds
end</code></pre>
</div>


<div class="markdown"><h2>Comparison</h2>
</div>


<div class="markdown"><h5>Statistics</h5>
<p>We compute here basic statistics, <em>i.e.</em> mean, median, min, max, standard deviation. </p>
<pre><code class="language-julia">statsPB &#61; basic_stat&#40;time_PB&#41;
statsPE &#61; basic_stat&#40;time_PE&#41;
statsLG &#61; basic_stat&#40;time_LG&#41;

factorPE &#61; mean&#40;time_PE&#41; / mean&#40;time_PB&#41;
factorLG &#61; mean&#40;time_LG&#41; / mean&#40;time_PB&#41;

# Write overall timings:
df &#61; DataFrame&#40;
	&#91;getfield&#40;j, i&#41; for i in fieldnames&#40;StatResults&#41;, j in &#91;statsPB, statsPE, statsLG&#93;&#93;,
	&#91;&quot;PlantBiophysics&quot;, &quot;plantecophys&quot;, &quot;LeafGasExchange&quot;&#93;
&#41;
insertcols&#33;&#40;df, 1, :Stat &#61;&gt; &#91;fieldnames&#40;StatResults&#41;...&#93;&#41;
CSV.write&#40;&quot;benchmark.csv&quot;, df&#41;

# Write timing for each sample:
CSV.write&#40;&quot;benchmark_full.csv&quot;,
	DataFrame&#40;
		&quot;package&quot; &#61;&gt; vcat&#40;
			&#91;
				repeat&#40;&#91;i.first&#93;, length&#40;i.second&#41;&#41; for i in &#91;
					&quot;PlantBiophysics&quot; &#61;&gt; time_PB,
					&quot;plantecophys&quot; &#61;&gt; time_PE,
					&quot;LeafGasExchange&quot; &#61;&gt; time_LG
				&#93;
			&#93;...
		&#41;,
		&quot;sample_time&quot; &#61;&gt; vcat&#40;time_PB, time_PE, time_LG&#41;
	&#41;
&#41;</code></pre>
</div>

<pre class='language-julia'><code class='language-julia'>df_res = CSV.read("benchmark.csv", DataFrame)</code></pre>
<table>
<tr>
<th></th>
<th>Stat</th>
<th>PlantBiophysics</th>
<th>plantecophys</th>
<th>LeafGasExchange</th>
</tr>
<tr>
<td>1</td>
<td>"mean"</td>
<td>1.21661e-6</td>
<td>0.0470207</td>
<td>0.0226825</td>
</tr>
<tr>
<td>2</td>
<td>"median"</td>
<td>1.083e-6</td>
<td>0.04569</td>
<td>0.0214673</td>
</tr>
<tr>
<td>3</td>
<td>"stddev"</td>
<td>1.20394e-6</td>
<td>0.00808695</td>
<td>0.0047238</td>
</tr>
<tr>
<td>4</td>
<td>"min"</td>
<td>5.41e-7</td>
<td>0.0346401</td>
<td>0.0192761</td>
</tr>
<tr>
<td>5</td>
<td>"max"</td>
<td>2.5e-5</td>
<td>0.172202</td>
<td>0.12373</td>
</tr>
</table>



<div class="markdown"><h5>Histogram plotting</h5>
</div>

<pre class='language-julia'><code class='language-julia'>Markdown.parse("""
We here display the computational time histogram of each package on the same scale in order to compare them: `PlantBiophysics.jl` (a), `plantecophys` (b) and `LeafGasExchange.jl` (c). The y-axis represents the density (_i.e._ reaching 0.3 means that 30% of the computed times are in this bar). Orange zone represents the interval [mean - standard deviation; mean + standard deviation]. Red dashed line represents the mean. Note that the x-axis is logarithmic.

```julia
fig = plot_benchmark_Makie(statsPB, statsPE, statsLG, time_PB, time_PE, time_LG)
    save("benchmark_each_time_steps.png", fig, px_per_unit=3)
```

![](https://github.com/VEZY/PlantBiophysics-paper/blob/main/tutorials/out/benchmark_each_time_steps.png?raw=true)

!!! note
    PlantBiophysics.jl is $(Int(round(df_res[1,:plantecophys] / df_res[1,:PlantBiophysics]))) times faster than plantecophys, and  $(Int(round(df_res[1,:LeafGasExchange] / df_res[1,:PlantBiophysics]))) times faster than LeafGasExchange.jl.

!!! warning
    This is the plot from the latest commit on &lt;https://github.com/VEZY/PlantBiophysics-paper/&gt;. If you want to make your own benchmarking, run the script that was used to perform it, but careful, it takes a long time to perform!
""")</code></pre>
<div class="markdown"><p>We here display the computational time histogram of each package on the same scale in order to compare them: <code>PlantBiophysics.jl</code> &#40;a&#41;, <code>plantecophys</code> &#40;b&#41; and <code>LeafGasExchange.jl</code> &#40;c&#41;. The y-axis represents the density &#40;<em>i.e.</em> reaching 0.3 means that 30&#37; of the computed times are in this bar&#41;. Orange zone represents the interval &#91;mean - standard deviation; mean &#43; standard deviation&#93;. Red dashed line represents the mean. Note that the x-axis is logarithmic.</p>
<pre><code class="language-julia">fig &#61; plot_benchmark_Makie&#40;statsPB, statsPE, statsLG, time_PB, time_PE, time_LG&#41;
    save&#40;&quot;benchmark_each_time_steps.png&quot;, fig, px_per_unit&#61;3&#41;</code></pre>
<p><img src="https://github.com/VEZY/PlantBiophysics-paper/blob/main/tutorials/out/benchmark_each_time_steps.png?raw&#61;true" alt="" /></p>
<div class="admonition is-note">
  <header class="admonition-header">
    Note
  </header>
  <div class="admonition-body">
    <p>PlantBiophysics.jl is 38649 times faster than plantecophys, and  18644 times faster than LeafGasExchange.jl.</p>
  </div>
</div>
<div class="admonition is-warning">
  <header class="admonition-header">
    Warning
  </header>
  <div class="admonition-body">
    <p>This is the plot from the latest commit on <a href="https://github.com/VEZY/PlantBiophysics-paper/">https://github.com/VEZY/PlantBiophysics-paper/</a>. If you want to make your own benchmarking, run the script that was used to perform it, but careful, it takes a long time to perform&#33;</p>
  </div>
</div>
</div>


<div class="markdown"><h1>References</h1>
</div>

<pre class='language-julia'><code class='language-julia'>"""
    StatResults(
        mean::AbstractFloat
    	median::AbstractFloat
    	stddev::AbstractFloat
    	min::AbstractFloat
    	max::AbstractFloat
    )

Structure to hold basic statistics of model performance.
"""
struct StatResults
    mean::AbstractFloat
    median::AbstractFloat
    stddev::AbstractFloat
    min::AbstractFloat
    max::AbstractFloat
end</code></pre>


<pre class='language-julia'><code class='language-julia'>begin
    function Base.show(io::IO, ::MIME"text/plain", m::StatResults)
        print(
            io,
            "Benchmark:",
            "\nMean time -&gt; ",
            m.mean,
            " ± ",
            m.stddev,
            "\nMedian time -&gt; ",
            m.median,
            "\nMinimum time -&gt; ",
            m.min,
            "\nMaximum time -&gt; ",
            m.max,
        )
    end
    Base.show(io::IO, m::StatResults) = print(io, m.mean, "(±", m.stddev, ')')


    md"""
    **Base.show**

    	Base.show(io::IO, m::StatResults)
    	Base.show(io::IO, ::MIME"text/plain", m::StatResults)

    Add a show method for our `StatResults` type.
    """
end</code></pre>
<div class="markdown"><p><strong>Base.show</strong></p>
<pre><code>Base.show&#40;io::IO, m::StatResults&#41;
Base.show&#40;io::IO, ::MIME&quot;text/plain&quot;, m::StatResults&#41;</code></pre>
<p>Add a show method for our <code>StatResults</code> type.</p>
</div>

<pre class='language-julia'><code class='language-julia'>"""
    basic_stat(df)

Compute basic statistics from the benchmarking
"""
function basic_stat(df)
    m = mean(df)
    med = median(df)
    std = Statistics.std(df)
    min = findmin(df)[1]
    max = findmax(df)[1]
    return StatResults(m, med, std, min, max)
end</code></pre>


<pre class='language-julia'><code class='language-julia'>function plot_benchmark_Makie(statsPB, statsPE, statsLG, time_PB, time_PE, time_LG)
    size_inches = (6.7, 5)
    size_pt = 72 .* size_inches
    bins = 220
    noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
    fig = Figure(
        backgroundcolor=RGBf(1, 1, 1),
        resolution=size_pt,
        font=noto_sans,
        fontsize=10,
    )
    ep = 1e-9
    extr = extrema(vcat(time_PB, time_PE, time_LG))
    interval = (extr[1] * 1e-1, extr[2])

    ga = fig[1, 1] = GridLayout()

    axa = Axis(
        ga[1, 1],
        title="(a) PlantBiophysics.jl",
        xscale=log10,
        titlealign=:left,
        titlesize=10,
    )
    stddevi = poly!(
        axa,
        Rect(max(ep, statsPB.mean - statsPB.stddev), 0.0, 2 * statsPB.stddev, 1),
        color=(:orange, 0.3),
        yautolimits=false,
    )
    moy = vlines!(axa, statsPB.mean; color=:red, linewidth=3, linestyle=:dot)
    hist!(axa, time_PB, normalization=:probability, bins=bins)
    # h = axa.finallimits[].widths[2]
    axislegend(
        axa,
        [stddevi, moy],
        ["95% confidence interval", "Mean"],
        "",
        position=:rb,
        orientation=:vertical,
        labelsize=8,
        framevisible=false,
    )
    xlims!(axa, interval)

    axb = Axis(
        ga[2, 1],
        title="(b) plantecophys",
        xscale=log10,
        ylabel="Density",
        titlealign=:left,
        titlesize=10,
    )
    stddevi = poly!(
        axb,
        Rect(statsPE.mean - statsPE.stddev, 0.0, 2 * statsPE.stddev, 1),
        color=(:orange, 0.3),
        yautolimits=false,
    )
    vlines!(axb, statsPE.mean; color=:red, linewidth=3, linestyle=:dot)
    hist!(axb, time_PE, normalization=:probability, bins=bins)
    xlims!(axb, interval)

    # axc = Axis(ga[3, 1], title="(c) LeafGasExchange.jl", yminorticks=IntervalsBetween(10),
    #     xscale=log10, xminorticks=IntervalsBetween(10), yminorgridvisible=true, yminorticksvisible=true,
    #     xminorgridvisible=true, xminorticksvisible=true, xlabel="Time (s)")
    axc = Axis(
        ga[3, 1],
        title="(c) LeafGasExchange.jl",
        xscale=log10,
        xlabel="Time (s)",
        titlealign=:left,
        titlesize=10,
    )
    stddevi = poly!(
        axc,
        Rect(statsLG.mean - statsLG.stddev, 0.0, 2 * statsLG.stddev, 1),
        color=(:orange, 0.3),
        yautolimits=false,
    )
    vlines!(axc, statsLG.mean; color=:red, linewidth=3, linestyle=:dot)
    hist!(axc, time_LG, normalization=:probability, bins=bins)
    xlims!(axc, interval)

    rowgap!(ga, 7)
    hidexdecorations!(axa, grid=false)
    hidexdecorations!(axb, grid=false)
    fig
end</code></pre>
<pre id='var-plot_benchmark_Makie' class='code-output documenter-example-output'>plot_benchmark_Makie (generic function with 1 method)</pre>
<div class='manifest-versions'>
<p>Built with Julia 1.8.5 and</p>
CSV 0.10.4<br>
DataFrames 1.3.4<br>
PlantBiophysics 0.3.0
</div>

<!-- PlutoStaticHTML.End -->
~~~