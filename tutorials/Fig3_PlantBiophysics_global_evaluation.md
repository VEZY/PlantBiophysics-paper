~~~
<!-- PlutoStaticHTML.Begin -->
<!--
    # This information is used for caching.
    [PlutoStaticHTML.State]
    input_sha = "b8a856f229984983d305b3027f45d87c0583af5f7d4f981489368b9a1dc1a3bd"
    julia_version = "1.8.5"
-->

<div class="markdown"><h1>Global evaluation of PlantBiophysics.jl: observations vs simulations</h1>
<p>This Pluto notebook presents the computation of Fig. 3 from the scientific article. The notebook does not compute anything because it would imply a dependency on R &#40;<a href="https://github.com/RemkoDuursma/plantecophys">plantecophys</a>&#41; and Python &#40;<a href="https://github.com/cropbox/LeafGasExchange.jl">LeafGasExchange.jl</a>&#41;, which is possible but not performant in Pluto an Github actions. Instead, we display the code and the results only. If you want to reproduce the results, execute the code provided here, or the code from the script provided <a href="https://github.com/VEZY/PlantBiophysics-paper/blob/main/tutorials/Fig3_PlantBiophysics_global_evaluation_noPluto.jl">here</a>.</p>
<h2>Importing the dependencies:</h2>
<h6>Note</h6>
<p>Make sure to have R installed on your computer first.</p>
</div>


<div class="markdown"><p>Loading the Julia packages:</p>
<pre><code class="language-julia">using CSV, Statistics, DataFrames, Downloads, Dates
using AlgebraOfGraphics, CairoMakie, Colors
using PlantBiophysics, PlantSimEngine, PlantMeteo, RCall, LeafGasExchange, Cropbox
using MonteCarloMeasurements</code></pre>
</div>

<pre class='language-julia'><code class='language-julia'>begin
    using CSV, Statistics, DataFrames, Downloads, Dates
    using AlgebraOfGraphics, CairoMakie, Colors
    using PlantBiophysics, PlantSimEngine, PlantMeteo
    using PlutoUI
end</code></pre>



<div class="markdown"><p>Loading the R package &#40;note the use of the <code>R&quot;&quot;&quot;</code> macro here that sends the code to R&#41;:</p>
<pre><code class="language-julia">R&quot;&quot;&quot;
library&#40;plantecophys&#41;
&quot;&quot;&quot;</code></pre>
</div>


<div class="markdown"><h2>Reading the data</h2>
</div>


<div class="markdown"><p>The data comes from Medlyn et al. &#40;2015&#41;, see <a href="https://figshare.com/articles/dataset/Tumbarumba_Gas_Exchange/1538079?file&#61;3402641">here</a> for more details.</p>
</div>

<pre class='language-julia'><code class='language-julia'>df = let 
    df_ = read_licor6400(Downloads.download("https://figshare.com/ndownloader/files/3402635"))

    # Computing DateTime from Date and Time:
    transform!(
        df_, 
        [:Date, :Time] =&gt; ((x, y) -&gt; Date.(x, dateformat"Y/m/d") .+ y) =&gt; :Date
    )

    # Initializing the columns:
    df_.VcMaxRef .= df_.JMaxRef .= df_.RdRef .= df_.TPURef .= df_.g0 .= df_.g1 .= df_.Tᵣ .= 0.0
    df_.AsimPB .= df_.EsimPB .= df_.TlsimPB .= df_.GssimPB .= 0.0
    df_.VcMaxRefPE .= df_.JMaxRefPE .= df_.RdRefPE .= df_.TPURefPE .= df_.g0PE .= df_.g1PE .= 0.0
    df_
end</code></pre>
<table>
<tr>
<th></th>
<th>Date</th>
<th>Time</th>
<th>Curve</th>
<th>Qflag</th>
<th>Site</th>
<th>Leaf Age</th>
<th>Chl a+b</th>
<th>Na</th>
<th>...</th>
</tr>
<tr>
<td>1</td>
<td>2001-11-14T09:38:00</td>
<td>09:38:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>2</td>
<td>2001-11-14T09:40:00</td>
<td>09:40:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>3</td>
<td>2001-11-14T09:42:00</td>
<td>09:42:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>4</td>
<td>2001-11-14T09:44:00</td>
<td>09:44:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>5</td>
<td>2001-11-14T09:48:00</td>
<td>09:48:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>6</td>
<td>2001-11-14T09:50:00</td>
<td>09:50:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>7</td>
<td>2001-11-14T09:52:00</td>
<td>09:52:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>8</td>
<td>2001-11-14T09:54:00</td>
<td>09:54:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>9</td>
<td>2001-11-14T09:56:00</td>
<td>09:56:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>10</td>
<td>2001-11-14T09:58:00</td>
<td>09:58:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>...</td>
</tr>
<tr>
<td>672</td>
<td>2002-05-10T15:12:00</td>
<td>15:12:00</td>
<td>48</td>
<td>1</td>
<td>5</td>
<td>1</td>
<td>"487.6722449"</td>
<td>2.20796</td>
<td></td>
</tr>
</table>



<div class="markdown"><h2>Fitting</h2>
<p>Photosynthesis and stomatal conductance parameters are fitted using PlantBiophysics.jl here, but it is also possible to use plantecophys instead &#40;see Julia script for more info&#41;. Comparison results are close with both.</p>
</div>

<pre class='language-julia'><code class='language-julia'>df_fit = let df_ = df
    for i in unique(df.Curve)
        dfi = filter(x -&gt; x.Curve == i, df_)
        sort!(dfi, :Cᵢ)
    
        g0, g1 = PlantSimEngine.fit(Medlyn, dfi)
        df_.g0[df.Curve.==i] .= g0
        df_.g1[df.Curve.==i] .= g1
    
        filter!(x -&gt; x.PPFD &gt; 1400.0, dfi)
    
        VcMaxRef, JMaxRef, RdRef, TPURef, Tᵣ = PlantSimEngine.fit(Fvcb, dfi)
        df_.VcMaxRef[df_.Curve.==i] .= VcMaxRef
        df_.JMaxRef[df_.Curve.==i] .= JMaxRef
        df_.RdRef[df_.Curve.==i] .= RdRef
        df_.TPURef[df_.Curve.==i] .= TPURef
        df_.Tᵣ[df_.Curve.==i] .= Tᵣ
    end
end</code></pre>



<div class="markdown"><h2>Simulation</h2>
<p>We set the wind velocity to 20 <span class="tex">$m.s^&#123;-1&#125;$</span> considering that the Licor-6400 chamber is well ventilated. The characteristic length is set to the square root of the chamber area in <span class="tex">$m^2$</span>. The leaf absorbtance and emissivity are set to the default value from plantecophys.</p>
</div>


<div class="markdown"><h3>PlantBiophysics</h3>
</div>

<pre class='language-julia'><code class='language-julia'>begin
    d = sqrt(df.Area[1]) / 100 # Characteristic dimension
    Wind = 20.0 # Wind, in m/s
    Leaf_abs = 0.86 # default from plantecophys
    emissivity = 0.95 # default from plantecophys
end</code></pre>
<pre id='var-emissivity' class='code-output documenter-example-output'>0.95</pre>

<pre class='language-julia'><code class='language-julia'>df_PB = let df = df  
    constants = Constants()
    atm_cols = keys(Atmosphere(T=25.0, Rh=0.5, Wind=10.0))
    for i in unique(df.Curve)
        dfi = filter(x -&gt; x.Curve == i, df)
    
        cols = fieldnames(Atmosphere)
        dfiMeteo = select(dfi, names(dfi, x -&gt; Symbol(x) in atm_cols))
        dfiMeteo.Wind .= Wind
        # Note that as we only use A-Ci curves, there is no NIR in the Licor6400
        dfiMeteo.Ri_SW_f .= dfi.PPFD .* Leaf_abs ./ (4.57)
        dfiMeteo.check .= false # Remove checks from Atmopshere (P &lt; 87kPa)
        meteo = Weather(dfiMeteo)
    
        leaf = ModelList(
            energy_balance=Monteith(
                aₛₕ=2,
                aₛᵥ=1,
                ε=emissivity, # Matching the value in plantecophys (https://github.com/RemkoDuursma/plantecophys/blob/c9749828041f10ca47c6691436678e0a5632cfb8/R/LeafEnergyBalance.R#L112)
                maxiter=100,
            ),
            photosynthesis=Fvcb(
                Tᵣ=dfi.Tᵣ[1],
                VcMaxRef=dfi.VcMaxRef[1],
                JMaxRef=dfi.JMaxRef[1],
                RdRef=dfi.RdRef[1],
                TPURef=dfi.TPURef[1],
            ),
            stomatal_conductance=Medlyn(dfi.g0[1], dfi.g1[1]),
            status=(Rₛ=meteo[:Ri_SW_f], sky_fraction=1.0, PPFD=dfi.PPFD, d=d)
        )
    
        run!(leaf, meteo)
        df.AsimPB[df.Curve.==i, :] = DataFrame(leaf).A
        df.EsimPB[df.Curve.==i, :] = DataFrame(leaf).λE ./ (meteo[:λ] * constants.Mₕ₂ₒ) * 1000
        df.TlsimPB[df.Curve.==i, :] = DataFrame(leaf).Tₗ
        df.GssimPB[df.Curve.==i, :] = DataFrame(leaf).Gₛ
    end
    df
end</code></pre>
<table>
<tr>
<th></th>
<th>Date</th>
<th>Time</th>
<th>Curve</th>
<th>Qflag</th>
<th>Site</th>
<th>Leaf Age</th>
<th>Chl a+b</th>
<th>Na</th>
<th>...</th>
</tr>
<tr>
<td>1</td>
<td>2001-11-14T09:38:00</td>
<td>09:38:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>2</td>
<td>2001-11-14T09:40:00</td>
<td>09:40:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>3</td>
<td>2001-11-14T09:42:00</td>
<td>09:42:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>4</td>
<td>2001-11-14T09:44:00</td>
<td>09:44:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>5</td>
<td>2001-11-14T09:48:00</td>
<td>09:48:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>6</td>
<td>2001-11-14T09:50:00</td>
<td>09:50:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>7</td>
<td>2001-11-14T09:52:00</td>
<td>09:52:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>8</td>
<td>2001-11-14T09:54:00</td>
<td>09:54:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>9</td>
<td>2001-11-14T09:56:00</td>
<td>09:56:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>10</td>
<td>2001-11-14T09:58:00</td>
<td>09:58:00</td>
<td>1</td>
<td>1</td>
<td>2</td>
<td>0</td>
<td>"434.0988776"</td>
<td>2.42908</td>
<td></td>
</tr>
<tr>
<td>...</td>
</tr>
<tr>
<td>672</td>
<td>2002-05-10T15:12:00</td>
<td>15:12:00</td>
<td>48</td>
<td>1</td>
<td>5</td>
<td>1</td>
<td>"487.6722449"</td>
<td>2.20796</td>
<td></td>
</tr>
</table>



<div class="markdown"><h3>plantecophys</h3>
<p>Below is presented the code for running the exact same simulation with the two other packages.</p>
<div class="admonition is-note">
  <header class="admonition-header">
    Note
  </header>
  <div class="admonition-body">
    <p>The code shown here for plantecophys and LeafGasExchange.jl is not executed, just displayed. If you want to reproduce the computation, please use the script instead.</p>
  </div>
</div>
</div>


<div class="markdown"><pre><code class="language-julia">begin
	if &#33;simulated
				
		R&quot;&quot;&quot;
		A_sim &#61; c&#40;&#41;
		E_sim &#61; c&#40;&#41;
		Tl_sim &#61; c&#40;&#41;
		Gs_sim &#61; c&#40;&#41;
		failed &#61;c&#40;&#41;
		&quot;&quot;&quot;

		df.AsimPE .&#61; df.EsimPE .&#61; df.TlsimPE .&#61; df.GssimPE .&#61; df.PEfailed .&#61; 0.
		for i in unique&#40;df.Curve&#41;
		    dfi &#61; filter&#40;x-&gt;x.Curve &#61;&#61; i,df&#41;
		    dfi &#61; dfi&#91;:,3:end&#93;
		    @rput dfi
		    Ca &#61; dfi.Cₐ
		    @rput Ca
		    R&quot;&quot;&quot;
		    VcMaxRef &#61; dfi&#36;VcMaxRefPE&#91;1&#93;
		    JMaxRef &#61; dfi&#36;JMaxRefPE&#91;1&#93;
		    TPURef &#61; dfi&#36;TPURefPE&#91;1&#93;
		    RdRef &#61; dfi&#36;RdRefPE&#91;1&#93;
		    g0 &#61; dfi&#36;g0PE&#91;1&#93;
		    g1 &#61; dfi&#36;g1PE&#91;1&#93;
		    res &#61; PhotosynEB&#40;Tair &#61; dfi&#36;T,Wind &#61; Wind,VPD&#61;dfi&#36;VPD,
		             Wleaf &#61; d,Ca &#61; Ca,  StomatalRatio &#61; 1,
		             LeafAbs &#61; 1.,gsmodel &#61; &quot;BBOpti&quot;,g0 &#61; g0, g1 &#61; g1,
		             EaV &#61; 58550.0,EdVC &#61; 2e&#43;05, delsC &#61; 629.26,
		             EaJ &#61; 29680.0,EdVJ &#61; 2e&#43;05,delsJ &#61; 631.88,
		             alpha &#61; 0.24,theta &#61; 0.7, Jmax &#61; JMaxRef, 
		             Vcmax &#61; VcMaxRef, TPU &#61; TPURef,Rd &#61; RdRef,
		             RH &#61; dfi&#36;Rh*100,PPFD&#61;dfi&#36;PPFD,
		             Patm &#61; dfi&#36;P,gk&#61;0.,
		             Tcorrect &#61; FALSE&#41;
		    A_sim &#61; append&#40;A_sim,res&#36;ALEAF&#41;
		    Gs_sim &#61; append&#40;Gs_sim,res&#36;GS&#41;
		    failed &#61; append&#40;failed,res&#36;failed&#41;
		    Tl_sim &#61; append&#40;Tl_sim,res&#36;Tleaf&#41;
		    E_sim &#61; append&#40;E_sim,res&#36;ELEAF&#41;
		    &quot;&quot;&quot;
		end
		@rget A_sim
		
		@rget failed
		@rget Gs_sim
		@rget E_sim
		@rget Tl_sim
		
		df.AsimPE .&#61; A_sim
		df.EsimPE .&#61; E_sim
		df.TlsimPE .&#61; Tl_sim
		df.GssimPE .&#61; Gs_sim
		df.PEfailed .&#61; failed
	end
end</code></pre>
</div>


<div class="markdown"><h3>LeafGasExchange.jl</h3>
</div>


<div class="markdown"><pre><code class="language-julia">begin
	df.AsimLG .&#61; df.EsimLG .&#61; df.TlsimLG .&#61; df.GssimLG .&#61; 0.0
	for i in unique&#40;df.Curve&#41;
	    dfi &#61; filter&#40;x -&gt; x.Curve &#61;&#61; i, df&#41;
	    configs &#61; &#91;&#93;
	    for i in 1:length&#40;dfi.T&#41;
	        config &#61; :Weather &#61;&gt; &#40;
	            PFD&#61;dfi.PPFD&#91;i&#93;,
	            CO2&#61;dfi.Cₐ&#91;i&#93;,
	            RH&#61;dfi.Rh&#91;i&#93; * 100,
	            T_air&#61;dfi.T&#91;i&#93;,
	            wind&#61;Wind,
	            P_air&#61;dfi.P&#91;i&#93;,
	            g0&#61;dfi.g0&#91;i&#93;,
	            g1&#61;1.57 * dfi.g1&#91;i&#93;,
	            Vc25&#61;dfi.VcMaxRef&#91;i&#93;,
	            Jm25&#61;dfi.JMaxRef&#91;i&#93;,
	            Rd25&#61;dfi.RdRef&#91;i&#93;,
	            Tp25&#61;dfi.TPURef&#91;i&#93;,
	            Ear&#61;46.39,
	            Haj&#61;29.68,
	            w&#61;d,
	            d&#61;d,
	            EaVc&#61;58.55,
	            ϵ&#61;emissivity,
	            Dh&#61;21.5,
	            α_s&#61;1 - Leaf_abs
	        &#41;
	        push&#33;&#40;configs, config&#41;
	    end
	    res &#61; simulate&#40;ModelC3MD; configs&#61;configs, nounit&#61;true&#41;
	    df.AsimLG&#91;df.Curve.&#61;&#61;i, :&#93; &#61; res.A_net
	    df.EsimLG&#91;df.Curve.&#61;&#61;i, :&#93; &#61; res.E
	    df.TlsimLG&#91;df.Curve.&#61;&#61;i, :&#93; &#61; res.T
	    df.GssimLG&#91;df.Curve.&#61;&#61;i, :&#93; &#61; res.gsc
	end
end</code></pre>
</div>


<div class="markdown"><h2>Stacking the results</h2>
<p>Here we stack the results in a long-format DataFrame for computing the statistics:</p>
</div>

<pre class='language-julia'><code class='language-julia'>begin
    meas =
        stack(
            select(
                df,
                [:Date =&gt; :Date, :Cₐ =&gt; :Cₐ, :A =&gt; :A, :Trmmol =&gt; :E, :Tₗ =&gt; :Tl, :Gₛ =&gt; :Gs]
            ),
            [:A, :E, :Tl, :Gs],
            [:Date, :Cₐ],
            value_name=:measured
        )
    
    sim_PB =
        stack(
            select(
                df_PB,
                [:Date =&gt; :Date, :AsimPB =&gt; :A, :EsimPB =&gt; :E, :TlsimPB =&gt; :Tl, :GssimPB =&gt; :Gs]
            ),
            [:A, :E, :Tl, :Gs],
            :Date,
            value_name=:simulated
        )
    sim_PB.origin .= "PlantBiophysics.jl"
    sim_PB
end</code></pre>
<table>
<tr>
<th></th>
<th>Date</th>
<th>variable</th>
<th>simulated</th>
<th>origin</th>
</tr>
<tr>
<td>1</td>
<td>2001-11-14T09:38:00</td>
<td>"A"</td>
<td>9.54724</td>
<td>"PlantBiophysics.jl"</td>
</tr>
<tr>
<td>2</td>
<td>2001-11-14T09:40:00</td>
<td>"A"</td>
<td>3.64789</td>
<td>"PlantBiophysics.jl"</td>
</tr>
<tr>
<td>3</td>
<td>2001-11-14T09:42:00</td>
<td>"A"</td>
<td>0.444311</td>
<td>"PlantBiophysics.jl"</td>
</tr>
<tr>
<td>4</td>
<td>2001-11-14T09:44:00</td>
<td>"A"</td>
<td>0.0</td>
<td>"PlantBiophysics.jl"</td>
</tr>
<tr>
<td>5</td>
<td>2001-11-14T09:48:00</td>
<td>"A"</td>
<td>19.2821</td>
<td>"PlantBiophysics.jl"</td>
</tr>
<tr>
<td>6</td>
<td>2001-11-14T09:50:00</td>
<td>"A"</td>
<td>30.654</td>
<td>"PlantBiophysics.jl"</td>
</tr>
<tr>
<td>7</td>
<td>2001-11-14T09:52:00</td>
<td>"A"</td>
<td>32.235</td>
<td>"PlantBiophysics.jl"</td>
</tr>
<tr>
<td>8</td>
<td>2001-11-14T09:54:00</td>
<td>"A"</td>
<td>33.4133</td>
<td>"PlantBiophysics.jl"</td>
</tr>
<tr>
<td>9</td>
<td>2001-11-14T09:56:00</td>
<td>"A"</td>
<td>29.3508</td>
<td>"PlantBiophysics.jl"</td>
</tr>
<tr>
<td>10</td>
<td>2001-11-14T09:58:00</td>
<td>"A"</td>
<td>20.3233</td>
<td>"PlantBiophysics.jl"</td>
</tr>
<tr>
<td>...</td>
</tr>
<tr>
<td>2688</td>
<td>2002-05-10T15:12:00</td>
<td>"Gs"</td>
<td>0.0819242</td>
<td>"PlantBiophysics.jl"</td>
</tr>
</table>



<div class="markdown"><p>Same for the other packages, though the data is not available in the notebook:</p>
</div>


<div class="markdown"><pre><code class="language-julia">begin 
	sim_LG &#61;
		    stack&#40;
		        select&#40;
		            df,
		            &#91;:Date &#61;&gt; :Date, :AsimLG &#61;&gt; :A, :EsimLG &#61;&gt; :E, :TlsimLG &#61;&gt; :Tl, :GssimLG &#61;&gt; :Gs&#93;
		        &#41;,
		        &#91;:A, :E, :Tl, :Gs&#93;,
		        :Date,
		        value_name&#61;:simulated
		    &#41;
		sim_LG.origin .&#61; &quot;LeafGasExchange.jl&quot;
		
		sim_PE &#61; stack&#40;
		    select&#40;
		        df,
		        &#91;:Date &#61;&gt; :Date, :AsimPE &#61;&gt; :A, :EsimPE &#61;&gt; :E, :TlsimPE &#61;&gt; :Tl, :GssimPE &#61;&gt; :Gs&#93;
		    &#41;,
		    &#91;:A, :E, :Tl, :Gs&#93;,
		    :Date,
		    value_name&#61;:simulated
		&#41;
		sim_PE.origin .&#61; &quot;plantecophys&quot;
end</code></pre>
</div>


<div class="markdown"><p>And finally we can stack the results together:</p>
</div>


<div class="markdown"><pre><code class="language-julia">begin 
	df_all &#61; vcat&#40;sim_PB, sim_LG, sim_PE&#41;
	df_res &#61; leftjoin&#40;df_all, meas, on&#61;&#91;:Date, :variable&#93;&#41;
end</code></pre>
</div>


<div class="markdown"><p>Note that the real code executed here is the following, because neither <code>sim_LG</code> or <code>sim_PE</code> are available:</p>
</div>

<pre class='language-julia'><code class='language-julia'>begin
    df_all = sim_PB
    df_res = leftjoin(df_all, meas, on=[:Date, :variable])
    # Filtering out the results for very low Cₐ:
    filter!(x -&gt; x.Cₐ &gt; 150, df_res)
    filter!(x -&gt; x.Cₐ &gt; 150, df)
    df_res
end</code></pre>
<table>
<tr>
<th></th>
<th>Date</th>
<th>variable</th>
<th>simulated</th>
<th>origin</th>
<th>Cₐ</th>
<th>measured</th>
</tr>
<tr>
<td>1</td>
<td>2001-11-14T09:38:00</td>
<td>"A"</td>
<td>9.54724</td>
<td>"PlantBiophysics.jl"</td>
<td>188.35</td>
<td>10.5883</td>
</tr>
<tr>
<td>2</td>
<td>2001-11-14T09:48:00</td>
<td>"A"</td>
<td>19.2821</td>
<td>"PlantBiophysics.jl"</td>
<td>374.267</td>
<td>21.7167</td>
</tr>
<tr>
<td>3</td>
<td>2001-11-14T09:50:00</td>
<td>"A"</td>
<td>30.654</td>
<td>"PlantBiophysics.jl"</td>
<td>760.9</td>
<td>33.0167</td>
</tr>
<tr>
<td>4</td>
<td>2001-11-14T09:52:00</td>
<td>"A"</td>
<td>32.235</td>
<td>"PlantBiophysics.jl"</td>
<td>955.4</td>
<td>34.8167</td>
</tr>
<tr>
<td>5</td>
<td>2001-11-14T09:54:00</td>
<td>"A"</td>
<td>33.4133</td>
<td>"PlantBiophysics.jl"</td>
<td>1155.4</td>
<td>36.85</td>
</tr>
<tr>
<td>6</td>
<td>2001-11-14T09:56:00</td>
<td>"A"</td>
<td>29.3508</td>
<td>"PlantBiophysics.jl"</td>
<td>1157.62</td>
<td>32.15</td>
</tr>
<tr>
<td>7</td>
<td>2001-11-14T09:58:00</td>
<td>"A"</td>
<td>20.3233</td>
<td>"PlantBiophysics.jl"</td>
<td>1168.7</td>
<td>23.9833</td>
</tr>
<tr>
<td>8</td>
<td>2001-11-14T10:00:00</td>
<td>"A"</td>
<td>9.77645</td>
<td>"PlantBiophysics.jl"</td>
<td>1184.58</td>
<td>10.7583</td>
</tr>
<tr>
<td>9</td>
<td>2001-11-14T10:02:00</td>
<td>"A"</td>
<td>5.16737</td>
<td>"PlantBiophysics.jl"</td>
<td>1195.37</td>
<td>6.26167</td>
</tr>
<tr>
<td>10</td>
<td>2001-11-14T10:04:00</td>
<td>"A"</td>
<td>2.59479</td>
<td>"PlantBiophysics.jl"</td>
<td>1197.85</td>
<td>2.07667</td>
</tr>
<tr>
<td>...</td>
</tr>
<tr>
<td>2228</td>
<td>2002-05-10T15:12:00</td>
<td>"Gs"</td>
<td>0.0819242</td>
<td>"PlantBiophysics.jl"</td>
<td>1199.09</td>
<td>0.0471338</td>
</tr>
</table>



<div class="markdown"><h2>Statistics</h2>
<p>Computing the statistics for each model:</p>
</div>

<pre class='language-julia'><code class='language-julia'>stats =
    combine(
        groupby(df_res, [:variable, :origin], sort=true),
        [:measured, :simulated] =&gt; ((x, y) -&gt; RMSE(x, y)) =&gt; :RMSE,
        [:measured, :simulated] =&gt; ((x, y) -&gt; nRMSE(x, y)) =&gt; :nRMSE,
        [:measured, :simulated] =&gt; ((x, y) -&gt; Bias(x, y)) =&gt; :Bias,
        [:measured, :simulated] =&gt; ((x, y) -&gt; nBias(x, y)) =&gt; :nBias,
        [:measured, :simulated] =&gt; ((x, y) -&gt; EF(x, y)) =&gt; :EF
    )</code></pre>
<table>
<tr>
<th></th>
<th>variable</th>
<th>origin</th>
<th>RMSE</th>
<th>nRMSE</th>
<th>Bias</th>
<th>nBias</th>
<th>EF</th>
</tr>
<tr>
<td>1</td>
<td>"A"</td>
<td>"PlantBiophysics.jl"</td>
<td>1.48</td>
<td>0.03</td>
<td>-0.0365</td>
<td>-0.0</td>
<td>0.98</td>
</tr>
<tr>
<td>2</td>
<td>"E"</td>
<td>"PlantBiophysics.jl"</td>
<td>0.48</td>
<td>0.05</td>
<td>0.2717</td>
<td>0.03</td>
<td>0.9</td>
</tr>
<tr>
<td>3</td>
<td>"Gs"</td>
<td>"PlantBiophysics.jl"</td>
<td>0.02</td>
<td>0.07</td>
<td>0.0074</td>
<td>0.02</td>
<td>0.88</td>
</tr>
<tr>
<td>4</td>
<td>"Tl"</td>
<td>"PlantBiophysics.jl"</td>
<td>0.52</td>
<td>0.02</td>
<td>0.3435</td>
<td>0.02</td>
<td>0.99</td>
</tr>
</table>



<div class="markdown"><h2>Plotting</h2>
</div>


<div class="markdown"><h3>Observations vs simulations &#40;Fig. 3&#41;</h3>
</div>

<pre class='language-julia'><code class='language-julia'>begin
    transparency_col = 0.8
    transparency_fill = 0.3
    color_pb = rgb(0, 0, 0, transparency_col)
    color_lg = rgb(223, 120, 97, transparency_col)
    color_pe = rgb(118, 84, 154, transparency_col)
    fill_pb = rgb(0, 0, 0, transparency_fill)
    fill_lg = rgb(223, 120, 97, transparency_fill)
    fill_pe = rgb(118, 84, 154, transparency_fill)
    stw = 1.5 # strokewidth
    ms = 7 # markersize
    noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
    legend_lab_size = 10
    size_inches = (10, 10)
    size_pt = 72 .* size_inches
    fig = Figure(
        font=noto_sans,
        resolution=size_pt,
        fontsize=12,
        xminorgridstyle=true
    )

    sideinfo1 = Label(fig[1:2, 1], "Simulations", rotation=pi / 2, fontsize=12)
    sideinfo2 = Label(fig[3, 2:3], "Observations", fontsize=12)

    # Assimilation
    axa = Axis(fig[1, 2], title="a) Net CO₂ assimilation (Aₙ)", aspect=1, titlealign=:left)
    xlims!(-10.0, 50.0)
    ylims!(-10.0, 50.0)

    ablines!(axa, 0, 1, color=(:grey, 0.4), linewidth=4)

    # LG = scatter!(
    #     axa, df.A, df.AsimLG,
    #     color=fill_lg,
    #     markersize=ms,
    #     strokecolor=color_lg,
    #     strokewidth=stw
    # )
    # PE = scatter!(axa, df.A, df.AsimPE,
    #     color=fill_pe,
    #     markersize=ms,
    #     strokecolor=color_pe,
    #     strokewidth=stw
    # )
    PB = scatter!(axa, df.A, df.AsimPB,
        color=fill_pb,
        markersize=ms,
        strokecolor=color_pb,
        strokewidth=stw
    )
    axislegend(
        axa,
        [PB],
        [
            "nRMSE: " * string(filter(x -&gt; x.variable == "A" && x.origin == "PlantBiophysics.jl", stats).nRMSE[1]),
        ],
        "",
        position=:rb,
        orientation=:vertical,
        labelsize=legend_lab_size,
        padding=0.0,
        framevisible=false
    )

    # Transpiration
    axb = Axis(fig[1, 3], title="b) Transpiration rate (E)", aspect=1, titlealign=:left)
    xlims!(-0.5, 10.0)
    ylims!(-0.5, 10.0)

    ablines!(axb, 0, 1, color=(:grey, 0.4), linewidth=4)

    # LG = scatter!(axb, df.Trmmol, df.EsimLG,
    #     color=fill_lg,
    #     markersize=ms,
    #     strokecolor=color_lg,
    #     strokewidth=stw
    # )
    # PE = scatter!(axb, df.Trmmol, df.EsimPE,
    #     color=fill_pe,
    #     markersize=ms,
    #     strokecolor=color_pe,
    #     strokewidth=stw
    # )
    PB = scatter!(axb, df.Trmmol, df.EsimPB,
        color=fill_pb,
        markersize=ms,
        strokecolor=color_pb,
        strokewidth=stw
    )
    axislegend(
        axb,
        [PB],
        [
            "nRMSE: " * string(filter(x -&gt; x.variable == "E" && x.origin == "PlantBiophysics.jl", stats).nRMSE[1]),
            # "nRMSE: " * string(filter(x -&gt; x.variable == "E" && x.origin == "plantecophys", stats).nRMSE[1]),
            # "nRMSE: " * string(filter(x -&gt; x.variable == "E" && x.origin == "LeafGasExchange.jl", stats).nRMSE[1]),
        ],
        "", position=:rb,
        orientation=:vertical,
        labelsize=legend_lab_size,
        padding=0.0,
        framevisible=false
    )

    # Stomatal conductance
    axc = Axis(
        fig[2, 2],
        title="(c) CO₂ stomatal conductance (Gₛ)",
        aspect=1, 
        titlealign=:left
    )
    xlims!(-0.05, 0.85)
    ylims!(-0.05, 0.85)

    ablines!(axc, 0, 1, color=(:grey, 0.4), linewidth=4)

    # LG = scatter!(axc, df.Gₛ, df.GssimLG,
    #     color=fill_lg,
    #     markersize=ms,
    #     strokecolor=color_lg,
    #     strokewidth=stw,
    #     label="LeafGasExchange.jl"
    # )
    # PE = scatter!(axc, df.Gₛ, df.GssimPE,
    #     color=fill_pe,
    #     markersize=ms,
    #     strokecolor=color_pe,
    #     label="plantecophys",
    #     strokewidth=stw,
    # )
    PB = scatter!(axc, df.Gₛ, df.GssimPB,
        color=fill_pb,
        markersize=ms,
        strokecolor=color_pb,
        strokewidth=stw,
        label="PlantBiophysics.jl"
    )

    axislegend(
        axc,
        [PB],
        [
            "nRMSE: " * string(filter(x -&gt; x.variable == "Gs" && x.origin == "PlantBiophysics.jl", stats).nRMSE[1]),
            # "nRMSE: " * string(filter(x -&gt; x.variable == "Gs" && x.origin == "plantecophys", stats).nRMSE[1]),
            # "nRMSE: " * string(filter(x -&gt; x.variable == "Gs" && x.origin == "LeafGasExchange.jl", stats).nRMSE[1]),
        ],
        "",
        position=:rb,
        orientation=:vertical,
        labelsize=legend_lab_size,
        padding=0.0,
        framevisible=false
    )


    # Leaf temperature
    axd = Axis(fig[2, 3], title="(d) Leaf temperature (Tₗ)", aspect=1, titlealign=:left)
    ablines!(axd, 0, 1, color=(:grey, 0.4), linewidth=4)

    # LG = scatter!(axd, df.Tₗ, df.TlsimLG,
    #     color=fill_lg,
    #     markersize=ms,
    #     strokecolor=color_lg,
    #     strokewidth=stw
    # )
    # PE = scatter!(axd, df.Tₗ, df.TlsimPE,
    #     color=fill_pe,
    #     markersize=ms,
    #     strokecolor=color_pe,
    #     strokewidth=stw
    # )
    PB = scatter!(axd, df.Tₗ, df.TlsimPB,
        color=fill_pb,
        markersize=ms,
        strokecolor=color_pb,
        strokewidth=stw
    )

    xlims!(10.0, 36.0)
    ylims!(10.0, 36.0)

    axislegend(
        axd, [PB],
        [
            "nRMSE: " * string(filter(x -&gt; x.variable == "Tl" && x.origin == "PlantBiophysics.jl", stats).nRMSE[1]),
            # "nRMSE: " * string(filter(x -&gt; x.variable == "Tl" && x.origin == "plantecophys", stats).nRMSE[1]),
            # "nRMSE: " * string(filter(x -&gt; x.variable == "Tl" && x.origin == "LeafGasExchange.jl", stats).nRMSE[1]),
        ],
        "", position=:rb,
        orientation=:vertical,
        labelsize=legend_lab_size,
        padding=0.0,
        framevisible=false
    )

    Legend(
        fig[4, 1:end],
        axc,
        orientation=:horizontal,
        framevisible=false,
        padding=0.0
    )
    fig
end</code></pre>
<img src="data:image/png;base64,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N27cKC4ulslk7u7ukG28KvT0sA8ZMqSkpERXYpCH/cyZM6mpqdeuXTOSbAN144Tj6dOn58+fv3r16tatW7du3bp9+/a2b30Zc+bM8fT0RAgRBMFgMPCsxhYWFrrfLWD8GAxGWFhYaGhoQEDAO++888cffzTutNVhvXv3TkxMbFwSHBycnp6Oq2ERQrNnz25oaBgyZEhISIiPj8/atWtxWvD06VPfRhITE/HwB51//etfcrl81qxZc+bM2b59+9OnT/EOTcrxzpGRkVOmTHn5y3nNNTQ0xMTEuLm5BQQE8Pn8urq6urq6P//809TUdN68eYaODrSVnh52R0fH5OTkxiVd/7BfvHgxNjbW09MTj0praGh4+et6WXQ3dezYseDgYPy5qKjI1NS07VubFRoaGhoa2upuNTU1+/btE4lEXC7XyspKKBSKRKJPP/203Regf2VlZUql0tBRtK6goMDQIbROqVSWlZV1yqmqq6sfPXokl8vbe2B5eblarW65vL6+vqampi1na+MvfLe0YcOGDRs2vGgrSZKhoaFOTk48Ho/P569YsWLjxo0jR47k8/kTJ058+PBhV4ZqKIWFhRRFGToK49XGB02vD3vbddnD3m37cFRUVFhbW+PPNjY2Uqm0oaEBz8TS6laEUGJiYpM/RouLi/GeLXypRCIpLCzs06fP119/ffTo0eLiYltb29mzZ69evbrlAw1CJpMxmUy1Wm3oQFohl8uN8O41oVKpZDJZ41+hDmOz2XZ2dlqttr1XbWJiolQqlUplC+UymYwkSRar9Qcfpgt7kYMHD37//fcsFovD4chksiNHjowfPx730Tt69KjuxQJAqywsLDrW+mZjY9OuciPRbRMOvEg0/qxSqQiCaPySbXkrQqi+vr64uLhxiUqloihKd9TzJBJJcXExTdNardbLy+uHH35wcHAQiUR4awsHGopWqyVJ0hi6LreMJEkjvHtNkCSJ76ehA2kF+Y9W96QoqgvieRUdO3YMITR+/HhTU9Pbt28/fPgwLi5u7NixmzdvhmwDgBZ024TD3t4+MjISfy4uLra2tm7851rLWxFC48ePHz9+fOOSr776CiHUo0ePZr+utra2pqZG12GYwWD069fPyDuOqdVq3PRj6EBaIZVKX3TbjQdOW40/ThaLRZJkW+KEJT9epLa2ls1mm5qaIoQGDx6cnZ3t7++v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">


<div class="markdown"><p><em>Figure 1. Measured &#40;x&#41; and simulated &#40;y&#41; net carbon assimilation &#40;a&#41;, net transpiration rate &#40;b&#41;, stomatal conductance for CO2 &#40;c&#41; and leaf temperature &#40;d&#41;. Data comes from Eucalyptus delegatensis &#40;Medlyn et al., 2015&#41;. All simulations were performed using a photosynthesis-stomatal conductance-energy balance coupled model with PlantBiophysics.jl. Grey line represents x&#61;y. All simulations were done with Ca &gt; 150 ppm</em></p>
</div>


<div class="markdown"><p>The full validation that includes all three packages would give the following plot:</p>
<p><img src="https://raw.githubusercontent.com/VEZY/PlantBiophysics-paper/main/tutorials/out/figure_global_simulation.png" alt="Full validation plot" /></p>
</div>


<div class="markdown"><h2>References</h2>
</div>

<pre class='language-julia'><code class='language-julia'>"""
    rgb(r, g, b, a)

Like `Colors.RGBA` but accepts colors in the 0-255 range.
"""
function rgb(r, g, b, a)
    return RGBA(r / 255, g / 255, b / 255, a)
end</code></pre>


<pre class='language-julia'><code class='language-julia'>"""
    RMSE(obs,sim)

Returns the Root Mean Squared Error between observations `obs` and simulations `sim`.
The closer to 0 the better.
"""
function RMSE(obs, sim, digits=2)
    return round(sqrt(sum((obs .- sim) .^ 2) / length(obs)), digits=digits)
end</code></pre>


<pre class='language-julia'><code class='language-julia'>"""
    nRMSE(obs,sim)

Returns the normalized Root Mean Squared Error between observations `obs` and simulations `sim`.
The closer to 0 the better.
"""
function nRMSE(obs, sim; digits=2)
    return round(sqrt(sum((obs .- sim) .^ 2) / length(obs)) / (findmax(obs)[1] - findmin(obs)[1]), digits=digits)
end</code></pre>


<pre class='language-julia'><code class='language-julia'>"""
    EF(obs,sim)

Returns the Efficiency Factor between observations `obs` and simulations `sim` using NSE (Nash-Sutcliffe efficiency) model.
More information can be found at https://en.wikipedia.org/wiki/Nash%E2%80%93Sutcliffe_model_efficiency_coefficient.
The closer to 1 the better.
"""
function EF(obs, sim, digits=2)
    SSres = sum((obs - sim) .^ 2)
    SStot = sum((obs .- mean(obs)) .^ 2)
    return round(1 - SSres / SStot, digits=digits)
end</code></pre>


<pre class='language-julia'><code class='language-julia'>"""
        Bias(obs,sim)

    Returns the bias between observations `obs` and simulations `sim`.
    The closer to 0 the better.
    """
function Bias(obs, sim, digits=4)
    return round(mean(sim .- obs), digits=digits)
end</code></pre>


<pre class='language-julia'><code class='language-julia'>"""
    nBias(obs,sim; digits = 2)

Returns the normalised bias (%) between observations `obs` and simulations `sim`.
The closer to 0 the better.
"""
function nBias(obs, sim; digits=2)
    return round(mean((sim .- obs)) / (findmax(obs)[1] - findmin(obs)[1]), digits=digits)
end</code></pre>


<pre class='language-julia'><code class='language-julia'>TableOfContents(title="📚 Table of Contents", indent=true, depth=4, aside=true)</code></pre>
<script>
	
const indent = true
const aside = true
const title_text = "📚 Table of Contents"
const include_definitions = false


const tocNode = html`<nav class="plutoui-toc">
	<header>
	 <span class="toc-toggle open-toc"></span>
	 <span class="toc-toggle closed-toc"></span>
	 ${title_text}
	</header>
	<section></section>
</nav>`

tocNode.classList.toggle("aside", aside)
tocNode.classList.toggle("indent", indent)


const getParentCell = el => el.closest("pluto-cell")

const getHeaders = () => {
	const depth = Math.max(1, Math.min(6, 4)) // should be in range 1:6
	const range = Array.from({length: depth}, (x, i) => i+1) // [1, ..., depth]
	
	const selector = [
		...(include_definitions ? [
			`pluto-notebook pluto-cell .pluto-docs-binding`, 
			`pluto-notebook pluto-cell assignee:not(:empty)`, 
		] : []),
		...range.map(i => `pluto-notebook pluto-cell h${i}`)
	].join(",")
	return Array.from(document.querySelectorAll(selector)).filter(el => 
		// exclude headers inside of a pluto-docs-binding block
		!(el.nodeName.startsWith("H") && el.closest(".pluto-docs-binding"))
	)
}


const document_click_handler = (event) => {
	const path = (event.path || event.composedPath())
	const toc = path.find(elem => elem?.classList?.contains?.("toc-toggle"))
	if (toc) {
		event.stopImmediatePropagation()
		toc.closest(".plutoui-toc").classList.toggle("hide")
	}
}

document.addEventListener("click", document_click_handler)


const header_to_index_entry_map = new Map()
const currently_highlighted_set = new Set()

const last_toc_element_click_time = { current: 0 }

const intersection_callback = (ixs) => {
	let on_top = ixs.filter(ix => ix.intersectionRatio > 0 && ix.intersectionRect.y < ix.rootBounds.height / 2)
	if(on_top.length > 0){
		currently_highlighted_set.forEach(a => a.classList.remove("in-view"))
		currently_highlighted_set.clear()
		on_top.slice(0,1).forEach(i => {
			let div = header_to_index_entry_map.get(i.target)
			div.classList.add("in-view")
			currently_highlighted_set.add(div)
			
			/// scroll into view
			/*
			const toc_height = tocNode.offsetHeight
			const div_pos = div.offsetTop
			const div_height = div.offsetHeight
			const current_scroll = tocNode.scrollTop
			const header_height = tocNode.querySelector("header").offsetHeight
			
			const scroll_to_top = div_pos - header_height
			const scroll_to_bottom = div_pos + div_height - toc_height
			
			// if we set a scrollTop, then the browser will stop any currently ongoing smoothscroll animation. So let's only do this if you are not currently in a smoothscroll.
			if(Date.now() - last_toc_element_click_time.current >= 2000)
				if(current_scroll < scroll_to_bottom){
					tocNode.scrollTop = scroll_to_bottom
				} else if(current_scroll > scroll_to_top){
					tocNode.scrollTop = scroll_to_top
				}
			*/
		})
	}
}
let intersection_observer_1 = new IntersectionObserver(intersection_callback, {
	root: null, // i.e. the viewport
  	threshold: 1,
	rootMargin: "-15px", // slightly smaller than the viewport
	// delay: 100,
})
let intersection_observer_2 = new IntersectionObserver(intersection_callback, {
	root: null, // i.e. the viewport
  	threshold: 1,
	rootMargin: "15px", // slightly larger than the viewport
	// delay: 100,
})

const render = (elements) => {
	header_to_index_entry_map.clear()
	currently_highlighted_set.clear()
	intersection_observer_1.disconnect()
	intersection_observer_2.disconnect()

		let last_level = `H1`
	return html`${elements.map(h => {
	const parent_cell = getParentCell(h)

		let [className, title_el] = h.matches(`.pluto-docs-binding`) ? ["pluto-docs-binding-el", h.firstElementChild] : [h.nodeName, h]

	const a = html`<a 
		class="${className}" 
		title="${title_el.innerText}"
		href="#${parent_cell.id}"
	>${title_el.innerHTML}</a>`
	/* a.onmouseover=()=>{
		parent_cell.firstElementChild.classList.add(
			'highlight-pluto-cell-shoulder'
		)
	}
	a.onmouseout=() => {
		parent_cell.firstElementChild.classList.remove(
			'highlight-pluto-cell-shoulder'
		)
	} */
		
		
	a.onclick=(e) => {
		e.preventDefault();
		last_toc_element_click_time.current = Date.now()
		h.scrollIntoView({
			behavior: 'smooth', 
			block: 'start'
		})
	}

	const row =  html`<div class="toc-row ${className} after-${last_level}">${a}</div>`
		intersection_observer_1.observe(title_el)
		intersection_observer_2.observe(title_el)
		header_to_index_entry_map.set(title_el, row)

	if(className.startsWith("H"))
		last_level = className
		
	return row
})}`
}

const invalidated = { current: false }

const updateCallback = () => {
	if (!invalidated.current) {
		tocNode.querySelector("section").replaceWith(
			html`<section>${render(getHeaders())}</section>`
		)
	}
}
updateCallback()
setTimeout(updateCallback, 100)
setTimeout(updateCallback, 1000)
setTimeout(updateCallback, 5000)

const notebook = document.querySelector("pluto-notebook")


// We have a mutationobserver for each cell:
const mut_observers = {
	current: [],
}

const createCellObservers = () => {
	mut_observers.current.forEach((o) => o.disconnect())
	mut_observers.current = Array.from(notebook.querySelectorAll("pluto-cell")).map(el => {
		const o = new MutationObserver(updateCallback)
		o.observe(el, {attributeFilter: ["class"]})
		return o
	})
}
createCellObservers()

// And one for the notebook's child list, which updates our cell observers:
const notebookObserver = new MutationObserver(() => {
	updateCallback()
	createCellObservers()
})
notebookObserver.observe(notebook, {childList: true})

// And finally, an observer for the document.body classList, to make sure that the toc also works when it is loaded during notebook initialization
const bodyClassObserver = new MutationObserver(updateCallback)
bodyClassObserver.observe(document.body, {attributeFilter: ["class"]})

// Hide/show the ToC when the screen gets small
let m = matchMedia("(max-width: 1000px)")
let match_listener = () => 
	tocNode.classList.toggle("hide", m.matches)
match_listener()
m.addListener(match_listener)

invalidation.then(() => {
	invalidated.current = true
	intersection_observer_1.disconnect()
	intersection_observer_2.disconnect()
	notebookObserver.disconnect()
	bodyClassObserver.disconnect()
	mut_observers.current.forEach((o) => o.disconnect())
	document.removeEventListener("click", document_click_handler)
	m.removeListener(match_listener)
})

return tocNode
</script>
<style>
@media not print {

.plutoui-toc {
	font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Oxygen-Sans, Cantarell, Helvetica, Arial, "Apple Color Emoji",
		"Segoe UI Emoji", "Segoe UI Symbol", system-ui, sans-serif;
	--main-bg-color: #fafafa;
	--pluto-output-color: hsl(0, 0%, 36%);
	--pluto-output-h-color: hsl(0, 0%, 21%);
	--sidebar-li-active-bg: rgb(235, 235, 235);
	--icon-filter: unset;
}

@media (prefers-color-scheme: dark) {
	.plutoui-toc {
		--main-bg-color: #303030;
		--pluto-output-color: hsl(0, 0%, 90%);
		--pluto-output-h-color: hsl(0, 0%, 97%);
		--sidebar-li-active-bg: rgb(82, 82, 82);
		--icon-filter: invert(1);
	}
}

.plutoui-toc.aside {
	color: var(--pluto-output-color);
	position: fixed;
	right: 1rem;
	top: 5rem;
	width: min(80vw, 300px);
	padding: 0.5rem;
	padding-top: 0em;
	/* border: 3px solid rgba(0, 0, 0, 0.15); */
	border-radius: 10px;
	/* box-shadow: 0 0 11px 0px #00000010; */
	max-height: calc(100vh - 5rem - 90px);
	overflow: auto;
	z-index: 40;
	background-color: var(--main-bg-color);
	transition: transform 300ms cubic-bezier(0.18, 0.89, 0.45, 1.12);
}

.plutoui-toc.aside.hide {
	transform: translateX(calc(100% - 28px));
}
.plutoui-toc.aside.hide section {
	display: none;
}
.plutoui-toc.aside.hide header {
	margin-bottom: 0em;
	padding-bottom: 0em;
	border-bottom: none;
}
}  /* End of Media print query */
.plutoui-toc.aside.hide .open-toc,
.plutoui-toc.aside:not(.hide) .closed-toc,
.plutoui-toc:not(.aside) .closed-toc {
	display: none;
}

@media (prefers-reduced-motion) {
  .plutoui-toc.aside {
	transition-duration: 0s;
  }
}

.toc-toggle {
	cursor: pointer;
    padding: 1em;
    margin: -1em;
    margin-right: -0.7em;
    line-height: 1em;
    display: flex;
}

.toc-toggle::before {
    content: "";
    display: inline-block;
    height: 1em;
    width: 1em;
    background-image: url("https://cdn.jsdelivr.net/gh/ionic-team/ionicons@5.5.1/src/svg/list-outline.svg");
	/* generated using https://dopiaza.org/tools/datauri/index.php */
    background-image: url("data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSI1MTIiIGhlaWdodD0iNTEyIiB2aWV3Qm94PSIwIDAgNTEyIDUxMiI+PHRpdGxlPmlvbmljb25zLXY1LW88L3RpdGxlPjxsaW5lIHgxPSIxNjAiIHkxPSIxNDQiIHgyPSI0NDgiIHkyPSIxNDQiIHN0eWxlPSJmaWxsOm5vbmU7c3Ryb2tlOiMwMDA7c3Ryb2tlLWxpbmVjYXA6cm91bmQ7c3Ryb2tlLWxpbmVqb2luOnJvdW5kO3N0cm9rZS13aWR0aDozMnB4Ii8+PGxpbmUgeDE9IjE2MCIgeTE9IjI1NiIgeDI9IjQ0OCIgeTI9IjI1NiIgc3R5bGU9ImZpbGw6bm9uZTtzdHJva2U6IzAwMDtzdHJva2UtbGluZWNhcDpyb3VuZDtzdHJva2UtbGluZWpvaW46cm91bmQ7c3Ryb2tlLXdpZHRoOjMycHgiLz48bGluZSB4MT0iMTYwIiB5MT0iMzY4IiB4Mj0iNDQ4IiB5Mj0iMzY4IiBzdHlsZT0iZmlsbDpub25lO3N0cm9rZTojMDAwO3N0cm9rZS1saW5lY2FwOnJvdW5kO3N0cm9rZS1saW5lam9pbjpyb3VuZDtzdHJva2Utd2lkdGg6MzJweCIvPjxjaXJjbGUgY3g9IjgwIiBjeT0iMTQ0IiByPSIxNiIgc3R5bGU9ImZpbGw6bm9uZTtzdHJva2U6IzAwMDtzdHJva2UtbGluZWNhcDpyb3VuZDtzdHJva2UtbGluZWpvaW46cm91bmQ7c3Ryb2tlLXdpZHRoOjMycHgiLz48Y2lyY2xlIGN4PSI4MCIgY3k9IjI1NiIgcj0iMTYiIHN0eWxlPSJmaWxsOm5vbmU7c3Ryb2tlOiMwMDA7c3Ryb2tlLWxpbmVjYXA6cm91bmQ7c3Ryb2tlLWxpbmVqb2luOnJvdW5kO3N0cm9rZS13aWR0aDozMnB4Ii8+PGNpcmNsZSBjeD0iODAiIGN5PSIzNjgiIHI9IjE2IiBzdHlsZT0iZmlsbDpub25lO3N0cm9rZTojMDAwO3N0cm9rZS1saW5lY2FwOnJvdW5kO3N0cm9rZS1saW5lam9pbjpyb3VuZDtzdHJva2Utd2lkdGg6MzJweCIvPjwvc3ZnPg==");
    background-size: 1em;
	filter: var(--icon-filter);
}

.aside .toc-toggle.open-toc:hover::before {
    background-image: url("https://cdn.jsdelivr.net/gh/ionic-team/ionicons@5.5.1/src/svg/arrow-forward-outline.svg");
	/* generated using https://dopiaza.org/tools/datauri/index.php */
    background-image: url("data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSI1MTIiIGhlaWdodD0iNTEyIiB2aWV3Qm94PSIwIDAgNTEyIDUxMiI+PHRpdGxlPmlvbmljb25zLXY1LWE8L3RpdGxlPjxwb2x5bGluZSBwb2ludHM9IjI2OCAxMTIgNDEyIDI1NiAyNjggNDAwIiBzdHlsZT0iZmlsbDpub25lO3N0cm9rZTojMDAwO3N0cm9rZS1saW5lY2FwOnJvdW5kO3N0cm9rZS1saW5lam9pbjpyb3VuZDtzdHJva2Utd2lkdGg6NDhweCIvPjxsaW5lIHgxPSIzOTIiIHkxPSIyNTYiIHgyPSIxMDAiIHkyPSIyNTYiIHN0eWxlPSJmaWxsOm5vbmU7c3Ryb2tlOiMwMDA7c3Ryb2tlLWxpbmVjYXA6cm91bmQ7c3Ryb2tlLWxpbmVqb2luOnJvdW5kO3N0cm9rZS13aWR0aDo0OHB4Ii8+PC9zdmc+");
}
.aside .toc-toggle.closed-toc:hover::before {
    background-image: url("https://cdn.jsdelivr.net/gh/ionic-team/ionicons@5.5.1/src/svg/arrow-back-outline.svg");
	/* generated using https://dopiaza.org/tools/datauri/index.php */
    background-image: url("data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSI1MTIiIGhlaWdodD0iNTEyIiB2aWV3Qm94PSIwIDAgNTEyIDUxMiI+PHRpdGxlPmlvbmljb25zLXY1LWE8L3RpdGxlPjxwb2x5bGluZSBwb2ludHM9IjI0NCA0MDAgMTAwIDI1NiAyNDQgMTEyIiBzdHlsZT0iZmlsbDpub25lO3N0cm9rZTojMDAwO3N0cm9rZS1saW5lY2FwOnJvdW5kO3N0cm9rZS1saW5lam9pbjpyb3VuZDtzdHJva2Utd2lkdGg6NDhweCIvPjxsaW5lIHgxPSIxMjAiIHkxPSIyNTYiIHgyPSI0MTIiIHkyPSIyNTYiIHN0eWxlPSJmaWxsOm5vbmU7c3Ryb2tlOiMwMDA7c3Ryb2tlLWxpbmVjYXA6cm91bmQ7c3Ryb2tlLWxpbmVqb2luOnJvdW5kO3N0cm9rZS13aWR0aDo0OHB4Ii8+PC9zdmc+");
}



.plutoui-toc header {
	display: flex;
	align-items: center;
	gap: .3em;
	font-size: 1.5em;
	/* margin-top: -0.1em; */
	margin-bottom: 0.4em;
	padding: 0.5rem;
	margin-left: 0;
	margin-right: 0;
	font-weight: bold;
	/* border-bottom: 2px solid rgba(0, 0, 0, 0.15); */
	position: sticky;
	top: 0px;
	background: var(--main-bg-color);
	z-index: 41;
}
.plutoui-toc.aside header {
	padding-left: 0;
	padding-right: 0;
}

.plutoui-toc section .toc-row {
	white-space: nowrap;
	overflow: hidden;
	text-overflow: ellipsis;
	padding: .1em;
	border-radius: .2em;
}

.plutoui-toc section .toc-row.H1 {
	margin-top: 1em;
}


.plutoui-toc.aside section .toc-row.in-view {
	background: var(--sidebar-li-active-bg);
}


	
.highlight-pluto-cell-shoulder {
	background: rgba(0, 0, 0, 0.05);
	background-clip: padding-box;
}

.plutoui-toc section a {
	text-decoration: none;
	font-weight: normal;
	color: var(--pluto-output-color);
}
.plutoui-toc section a:hover {
	color: var(--pluto-output-h-color);
}

.plutoui-toc.indent section a.H1 {
	font-weight: 700;
	line-height: 1em;
}

.plutoui-toc.indent section .after-H2 a { padding-left: 10px; }
.plutoui-toc.indent section .after-H3 a { padding-left: 20px; }
.plutoui-toc.indent section .after-H4 a { padding-left: 30px; }
.plutoui-toc.indent section .after-H5 a { padding-left: 40px; }
.plutoui-toc.indent section .after-H6 a { padding-left: 50px; }

.plutoui-toc.indent section a.H1 { padding-left: 0px; }
.plutoui-toc.indent section a.H2 { padding-left: 10px; }
.plutoui-toc.indent section a.H3 { padding-left: 20px; }
.plutoui-toc.indent section a.H4 { padding-left: 30px; }
.plutoui-toc.indent section a.H5 { padding-left: 40px; }
.plutoui-toc.indent section a.H6 { padding-left: 50px; }


.plutoui-toc.indent section a.pluto-docs-binding-el,
.plutoui-toc.indent section a.ASSIGNEE
	{
	font-family: JuliaMono, monospace;
	font-size: .8em;
	/* background: black; */
	font-weight: 700;
    font-style: italic;
	color: var(--cm-var-color); /* this is stealing a variable from Pluto, but it's fine if that doesnt work */
}
.plutoui-toc.indent section a.pluto-docs-binding-el::before,
.plutoui-toc.indent section a.ASSIGNEE::before
	{
	content: "> ";
	opacity: .3;
}
</style>

<div class='manifest-versions'>
<p>Built with Julia 1.8.5 and</p>
AlgebraOfGraphics 0.6.14<br>
CSV 0.10.9<br>
CairoMakie 0.10.2<br>
Colors 0.12.10<br>
DataFrames 1.5.0<br>
Downloads 1.6.0<br>
PlantBiophysics 0.9.0<br>
PlantMeteo 0.3.0<br>
PlantSimEngine 0.5.1<br>
PlutoUI 0.7.50
</div>

<!-- PlutoStaticHTML.End -->
~~~