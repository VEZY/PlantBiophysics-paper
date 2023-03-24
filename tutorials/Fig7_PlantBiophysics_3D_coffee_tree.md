~~~
<!-- PlutoStaticHTML.Begin -->
<!--
    # This information is used for caching.
    [PlutoStaticHTML.State]
    input_sha = "20d362bd18ef7ebc0ed4a9f806e6175b208df916323980fb6bca97b0e810a56b"
    julia_version = "1.8.5"
-->

<div class="markdown"><h1>PlantBiophysics.jl 3D global tree simulation</h1>
<p>This Pluto notebook presents the computation of Fig. 7 from the scientific article. It displays leaf temperature on a 3D coffee tree simulated by PlantBiophysics.jl. Non-Pluto Julia script is also available &#40;see <a href="https://github.com/VEZY/PlantBiophysics-paper/blob/main/tutorials/Fig7_PlantBiophysics_3D_coffee_tree_noPluto.jl">here</a>&#41;.</p>
<h2>Importing the dependencies:</h2>
<p>Loading the Julia packages:</p>
</div>

<pre class='language-julia'><code class='language-julia'>begin
    using PlantBiophysics, PlantGeom, PlantMeteo, PlantSimEngine
    using BenchmarkTools
    using Dates, DataFrames, CSV, Statistics
    using MultiScaleTreeGraph
    using PlutoUI
    nothing
end</code></pre>



<div class="markdown"><h2>Reading data</h2>
<h4>MTG file</h4>
</div>


<div class="markdown"><p>julia<code>mtg &#61; read_opf&#40;&quot;coffee.opf&quot;&#41;</code></p>
</div>


<div class="markdown"><h4>Meteorological data</h4>
</div>


<div class="markdown"><pre><code>weather &#61; PlantMeteo.read_weather&#40;&quot;meteo.csv&quot;,
    		:temperature &#61;&gt; :T,
    		:relativeHumidity &#61;&gt; &#40;x -&gt; x ./ 100&#41; &#61;&gt; :Rh,
    		:wind &#61;&gt; :Wind,
    		:atmosphereCO2_ppm &#61;&gt; :Cₐ,
    		date_format &#61; DateFormat&#40;&quot;yyyy/mm/dd&quot;&#41;
&#41;</code></pre>
</div>


<div class="markdown"><h4>Models list</h4>
</div>


<div class="markdown"><pre><code class="language-julia">file &#61; &quot;plant_coffee.yml&quot;
models &#61; read_model&#40;file&#41;

to_initialize&#40;models</code></pre>
</div>


<div class="markdown"><h2>Adding light interception data: incident radiation</h2>
</div>


<div class="markdown"><pre><code class="language-julia">transform&#33;&#40;
    mtg,
    &#91;:Ra_PAR_f, :Ra_NIR_f&#93; &#61;&gt; &#40;&#40;x, y&#41; -&gt; x &#43; y * 1.2&#41; &#61;&gt; :Rᵢ, # This would be the incident radiation
    &#91;:Ra_PAR_f, :Ra_NIR_f&#93; &#61;&gt; &#40;&#40;x, y&#41; -&gt; x &#43; y&#41; &#61;&gt; :Rₛ,
    :Ra_PAR_f &#61;&gt; &#40;x -&gt; x * 4.57&#41; &#61;&gt; :PPFD,
    &#40;x -&gt; 0.3&#41; &#61;&gt; :d,
    ignore_nothing &#61; true
&#41;</code></pre>
</div>


<div class="markdown"><h2>Running the simulation</h2>
</div>


<div class="markdown"><pre><code class="language-julia">run&#33;&#40;mtg, models, weather&#41;</code></pre>
</div>


<div class="markdown"><h2>Quickly benchmarking</h2>
<p>For sake of simplicity, we benchmark here in Sequential mode &#40;i.e. no parallelization&#41;.</p>
</div>


<div class="markdown"><pre><code class="language-julia">B &#61; @benchmark run&#33;&#40;&#36;mtg, &#36;models, &#36;weather,executor&#61;SequentialEx&#40;&#41;&#41;</code></pre>
</div>


<div class="markdown"><h2>Plotting the result in 3D</h2>
</div>


<div class="markdown"><pre><code class="language-julia">begin
	transform&#33;&#40;
    mtg,
    :Tₗ &#61;&gt; &#40;x -&gt; x&#91;1&#93;&#41; &#61;&gt; :Tₗ_1,
    ignore_nothing &#61; true
&#41;

	f, ax, p &#61; viz&#40;mtg, color &#61; :Tₗ_1&#41;
	colorbar&#40;f&#91;1, 2&#93;, p&#41;
	f
end</code></pre>
</div>


<div class="markdown"><p><img src="https://github.com/VEZY/PlantBiophysics-paper/blob/main/tutorials/out/3d_coffee.png?raw&#61;true" alt="" /></p>
</div>
<div class='manifest-versions'>
<p>Built with Julia 1.8.5 and</p>
BenchmarkTools 1.3.2<br>
CSV 0.10.9<br>
DataFrames 1.5.0<br>
MultiScaleTreeGraph 0.9.0<br>
PlantBiophysics 0.9.0<br>
PlantGeom 0.5.1<br>
PlantMeteo 0.3.0<br>
PlantSimEngine 0.5.1<br>
PlutoUI 0.7.50
</div>

<!-- PlutoStaticHTML.End -->
~~~