To Use the Madgraph Filters:
1. Place the python script in the PLUGIN folder in madrgaph
2. Change the name to user_filter.py
3. When generating a process add the flag --diagram_filter

Explanation of different filters
  user_filter.py 
    Used when generating gg H -> 4l processes
      Options in code to select for ZZ Zgam and GamGam vertices only
  user_filter_VH.py
    Used when generating pp -> H 2l processes
      Options in code to select for ZZ Zgam and GamGam vertices only
