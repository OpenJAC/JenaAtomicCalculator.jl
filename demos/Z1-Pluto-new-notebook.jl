### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ f18a9488-13b8-11f0-342b-f99afaf3ac19
begin
	using Pkg
	Pkg.develop(path="/home/fritzsch/fri/JAC.jl")
end

# ╔═╡ fd4ef417-311f-43aa-a601-bf8a9d66b890
using JAC

# ╔═╡ 79374db5-1fbd-4499-b8d9-c7b840bdee52
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(150px, 10%);
    	padding-right: max(150px, 10%);
	}
</style>
"""

# ╔═╡ 8d491602-2e68-4ce9-b949-859f91cfee9d
md"""
# Compute the low-lying levels of C``^{2+}`` 
"""

# ╔═╡ 75061d53-4fb7-45a6-b8bb-c199d52479da
md"""
Let us first invoke JAC, either from the own source-code basis or simply by **using JAC**
"""

# ╔═╡ 956d7bd6-d59f-4daa-88b6-9b8879918d79
md"""
The low-lying levels (level structure) of beryllium-like ions, and especially of C``^{2+}``, has been calculated in many case studies
in the literature. While the level structure of these ions is still quite simple, it exhibits a considerable admixture of the 
``2s^22p^2`` configuration already for the ``1s^{2}2s^{2}`` ``^{1}S_{0}`` ground level.

We here show how the low-lying levels of C``^{2+}`` can be readily calculated in JAC by either following the default settings or
by specifying further details for both, the SCF and configuration-interaction (CI) computations. As usual, we first need to 
specify a radial grid as well as the nuclear model for the subsequent computations:
"""

# ╔═╡ b0f8586c-db7e-4872-b4e3-ff6beafb4c6c


# ╔═╡ 4c2b4d2a-4b34-4f6a-9690-6efc9f2d3fea


# ╔═╡ 759da271-0adf-4ae7-9642-8620c4c34690


# ╔═╡ d70082f6-fa88-4b54-8306-46efe25138fc
md"""
# Still under construction !!!
"""

# ╔═╡ Cell order:
# ╟─79374db5-1fbd-4499-b8d9-c7b840bdee52
# ╠═8d491602-2e68-4ce9-b949-859f91cfee9d
# ╟─75061d53-4fb7-45a6-b8bb-c199d52479da
# ╟─f18a9488-13b8-11f0-342b-f99afaf3ac19
# ╠═fd4ef417-311f-43aa-a601-bf8a9d66b890
# ╠═956d7bd6-d59f-4daa-88b6-9b8879918d79
# ╠═b0f8586c-db7e-4872-b4e3-ff6beafb4c6c
# ╠═4c2b4d2a-4b34-4f6a-9690-6efc9f2d3fea
# ╠═759da271-0adf-4ae7-9642-8620c4c34690
# ╟─d70082f6-fa88-4b54-8306-46efe25138fc
