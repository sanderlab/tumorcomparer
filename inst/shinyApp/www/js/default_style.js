[
   {"selector":"node", "css": {
       "text-valign": "center",
       "text-halign": "center",
       "border-color": "black",
       "content": "data(label)",
       "background-color": "data(color)",
       "border-width": "0px",
       "width": "mapData(degree, 0, 10, 20, 40)",
       "height": "mapData(degree, 0, 10, 20, 40)"
       }},
  
    {"selector": "node:selected", "css": {
       "overlay-opacity": 0.8,
       "overlay-padding": 15,
       "overlay-color": "magenta"
    }},

    {"selector": "edge", "css": {
        "curve-style": "straight"
    }},

    {"selector": "edge[interactionType='pp']", "css": {
        "line-color": "grey",
        "target-arrow-shape": "triangle",
        "target-arrow-color": "grey",
        "arrow-scale": 2
    }},
    
    {"selector": "edge[interactionType='interacts_with']", "css": {
        "line-color": "grey",
        "target-arrow-shape": "triangle",
        "target-arrow-color": "grey",
        "arrow-scale": 2
    }},

    {"selector": "edge[interactionType='catalysis-precedes']", "css": {
        "line-color": "lightblue",
        "target-arrow-shape": "triangle",
        "target-arrow-color": "lightblue",
        "arrow-scale": 2
    }}
]
