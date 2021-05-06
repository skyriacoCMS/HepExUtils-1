import madgraph.core.drawing as drawing

def remove_diag(diag, model):
    """remove all diagram with quark in T-channel"""

    #diag is a Direct Accyclic Graph of the Feynman diagram
    #convert it to a full graph (the one we used for the plotting)
    #In that representation each vertex is associated to all the legs/particle and we can move quite freely inside the diagram
    # The structure representation is simpler and should help to apply complex filtering

    draw = drawing.FeynmanDiagram(diag, model)
    draw.load_diagram()
    # Diagram content  can be loop over three list:
    #  - all vertex of diagram (note that for each external particle we have a vertex with a single particle attached to it)
    #    for v in draw.vertexList:
    #  - all vertex corresponding to the initial/final state particles
    #    for v in draw.initial_vertex:
    # - all the particles (initial states / final states and propagator) of the diagram
    #    for p in draw.lineList
    #
    # All vertex can be assigned to a level by running
    draw.define_level()
    #      0: correspond to the initial state vertex
    #      1: correspond to the T-channel level
    #   you can use draw._debug_level() to text representation in that case
    #   For a single line the "begin" level will always be one level lower than the "end" level
    #    BUT for t-channel where both are equal and set to 1
    #print(draw.load_diagram())
    #
    # All vertex are of type VertexPoint
    #   They have only two relevant attributes
    #      self.lines : all the particle connecting to this vertex
    #      self.level : level of this vertex (see above)
    #
    # All particles are of type FeynmanLine
    #   They have the following relevant attributes
    #       self.begin: vertex associated to the beginning of the particles
    #       self.end:  vertex associated to the end of the particles
    #       self.get('id'): pdg code
    #       self.get('number'): number associated to the original DAG
    #       other attributes are ['polarization', 'number', 'onshell', 'state', 'model', 'from_group', 'loop_line']

    #remove all diagrams with 2 Higgs vertices
    # remove all 4 point interaction
    for v in draw.vertexList:
        if len(v.lines) > 3:
            return True
    ## Note you can only filter out one of these
    RemoveHZZ=True
    RemoveHZgam=True
    RemoveHgamgam=False
    ###
    #### extraneous particles and Hll vertices and H H vertices
    for v in draw.vertexList:
      hasL=False
      hasH=False
      for p in v.lines:
        ##HLL check
        if (abs(p.id) == 11 or abs(p.id) == 13) and hasH and hasL:
          return True
        ##Hqq check
        if (abs(p.id) == 1 or abs(p.id) == 2 or abs(p.id) == 3 or abs(p.id) == 4) and hasH:
          return True
        if (p.id == 25) and hasL:
          return True
        if p.id == 32:
          return True
        if abs(p.id) == 11 or abs(p.id) == 13:
          hasL=True
        if p.id == 25:
          hasH=True
        ### Remove New Higgs and Z prime
        if p.id == 9000005:
          return True
        if p.id == 9000008:
          return True

    print(draw._debug_level())
    for v in draw.vertexList:
      HasZ=False
      HasH=False
      HasGam=False
      for p in v.lines:
        if abs(p.id) == 22:
          HasGam=True
        if p.id == 25:
          HasH=True
        if abs(p.id) == 23:
          HasZ=True
      if HasZ and HasH and not HasGam:
        if RemoveHZZ:
          return True
      if HasZ and HasH and HasGam:
        if RemoveHZgam:
          return True
      if HasGam and HasH and not HasZ:
        if RemoveHgamgam:
          return True
    # check if initial/final  states quark are emitting a gluon
    #for external in draw.initial_vertex: # vertex with only one leg for initial/final state particles
        # external is VertexPoint object
        #p = external.lines[0] # FeynmanLine object
        #if abs(p.id) < 7: # PDG code less than 7 means a SM quark
        #    other_vertex = p.end # VertexPoint object
    return False
