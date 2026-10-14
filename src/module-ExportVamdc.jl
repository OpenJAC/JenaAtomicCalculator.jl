using Printf

"""
`module  JAC.ExportVamdc`
    ... a submodel of JAC that writes an `ExportDataSet.DataSet` as a MINIMAL, VAMDC/XSAMS-SHAPED XML
        sketch, not a schema-valid XSAMS document.

        SCOPE, DELIBERATELY AND SEVERELY NARROW. XSAMS (the VAMDC exchange schema) is a large, strict XML
        schema with its own controlled vocabularies, species identifiers (InChI-like atomic/ionic IDs) and
        deeply nested process structures; a schema-valid writer is its own substantial undertaking, far
        beyond what this module attempts. What follows is a hand-written XML sketch using XSAMS-like element
        names (Species, AtomicState, RadiativeTransition, ...) that shows the SHAPE of the mapping from
        JAC's own data to VAMDC's -- levels and transitions in a recognisable nested structure -- and
        nothing more. IT HAS NOT BEEN VALIDATED AGAINST THE REAL XSAMS XSD, carries no namespace
        declarations a real VAMDC node would require, and is not meant to be ingested by VAMDC tooling.
        See [[project_challenges_list]] S24 for what a genuine XSAMS writer would still need.
"""
module ExportVamdc

using Printf, ..ExportDataSet


"""
`ExportVamdc.writeMinimalXsams(stream::IO, data::ExportDataSet.DataSet)`
    ... to write `data` as the minimal, non-schema-valid XML sketch this module's own docstring describes.
        Nothing is returned.
"""
function  writeMinimalXsams(stream::IO, data::ExportDataSet.DataSet)
    println(stream, "<!-- NOT a schema-valid XSAMS document -- see ExportVamdc's own docstring. -->")
    println(stream, "<XSAMSData>")
    println(stream, "  <Species>")
    @printf(stream, "    <Atom element=\"%s\" ionCharge=\"%d\" nuclearCharge=\"%d\">\n",
                     data.element, data.ionCharge, data.nuclearCharge)
    println(stream, "      <AtomicStates>")
    for  lev in data.levels
        @printf(stream, "        <AtomicState stateID=\"S%d\">\n", lev.index)
        @printf(stream, "          <Label>%s</Label>\n", lev.label)
        @printf(stream, "          <TotalAngularMomentum>%3.1f</TotalAngularMomentum>\n", lev.J)
        @printf(stream, "          <Energy unit=\"1/cm\">%.3f</Energy>\n", lev.energy)
        println(stream, "        </AtomicState>")
    end
    println(stream, "      </AtomicStates>")
    println(stream, "      <RadiativeTransitions>")
    for  tr in data.transitions
        println(stream, "        <RadiativeTransition>")
        @printf(stream, "          <UpperState stateRef=\"S%d\"/>\n", tr.upperIndex)
        @printf(stream, "          <LowerState stateRef=\"S%d\"/>\n", tr.lowerIndex)
        @printf(stream, "          <TransitionProbabilityA unit=\"1/s\">%.4e</TransitionProbabilityA>\n", tr.aValue)
        println(stream, "        </RadiativeTransition>")
    end
    println(stream, "      </RadiativeTransitions>")
    println(stream, "    </Atom>")
    println(stream, "  </Species>")
    println(stream, "</XSAMSData>")

    return( nothing )
end

end # module
