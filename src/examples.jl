function example1()
    parent1 = Crystal("IRMOF-1.cif")
	infer_bonds!(parent1, true)
	query1 = moiety("2-!-p-phenylene.xyz")
	replacement1 = moiety("2-acetylamido-p-phenylene.xyz")
	search = query1 ∈ parent1
	example1 = substructure_replace(search, replacement1,
		nb_loc=Int(nb_locations(search)/4))
    view_crystal(example1)
end


function example2()
    parent = Crystal("IRMOF-1_noH.cif")
	infer_bonds!(parent, true)
	query = moiety("1,4-C-phenylene_noH.xyz")
	replacement = moiety("1,4-C-phenylene.xyz")
	repaired_xtal = replace(parent, query => replacement)
	view_crystal(repaired_xtal)
end


function example3()
    parent = Crystal("SIFSIX-2-Cu-i.cif", check_overlap=false)
	infer_bonds!(parent, true)
	q = moiety("disordered_ligand!.xyz")
	r = moiety("4-pyridyl.xyz")
	repaired = replace(parent, q => r)
	search = substructure_search(moiety("acetylene.xyz"), repaired, disconnected_component=true)
	active = substructure_replace(search, nothing)
	view_crystal(active)
end


function example4()
    parent = Crystal("UiO-66.cif")
	infer_bonds!(parent, true)
	query = moiety("BDC.xyz")
	replacement = moiety("formate_caps.xyz")
	with_defect = replace(parent, query => replacement; loc=[1, 3, 9])
	view_crystal(with_defect)
end
