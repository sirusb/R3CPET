<style>

.node {
  stroke: #fff;
  stroke-width: 1.5px;
}

.link {
  stroke: #999;
  stroke-opacity: .6;
}

</style>

<!--<script src="d3.min.js"></script>-->
<script type="text/javascript">var networkOutputBinding = new Shiny.OutputBinding();
  $.extend(networkOutputBinding, {
    find: function(scope) {
      return $(scope).find('.shiny-network-output');
    },
    renderValue: function(el, data) {
      
      //format nodes object
      var nodes = new Array();
      for (var i = 0; i < data.names.length; i++){
        nodes.push({"name": data.names[i]})
      }


      var width = 800;
      var height = 600;
    
      var lin = data.links
      var force = d3.layout.force()
        .nodes(nodes)
        .links(lin)
        .charge(-3000)
        .friction(0.6)
        .gravity(0.6)
        .linkDistance(40)        
        .size([width, height])
        .start();
      
      //remove the old graph
      var svg = d3.select(el).select("svg");      
      svg.remove();
      
      $(el).html("");
      
      //append a new one
      svg = d3.select(el).append("svg");
      
      svg.attr("width", width)
        .attr("height", height);
    
      var link = svg.selectAll("line.link")
          .data(lin)
        .enter().append("line")
          .attr("class", "link")
          .style("stroke-width", function(d) { return Math.sqrt(d.value); });
    
      var node = svg.selectAll("circle.node")
          .data(nodes)
        .enter().append("circle")
          .attr("class", "node")
          .attr("r", 10)
          .style("fill", function(d) {return "#feb24c"})
          //.style("fill", function(d) { return color(d.group); })
          .call(force.drag);
      node.append("title")
        .text(function(d) { return d.name; });
    
    var texts = svg.selectAll("text.label")
                .data(nodes)
                .enter().append("text")
                .attr("class", "label")
                .attr("fill", "black")
                .text(function(d) {  return d.name;  });
        
      force.on("tick", function() {
        link.attr("x1", function(d) { return d.source.x; })
            .attr("y1", function(d) { return d.source.y; })
            .attr("x2", function(d) { return d.target.x; })
            .attr("y2", function(d) { return d.target.y; });
    
        node.attr("cx", function(d) { return d.x; })
            .attr("cy", function(d) { return d.y; });
          
        texts.attr("transform", function(d) {
            return "translate(" + d.x + "," + d.y + ")";
        });  
      });
      
    }
  });
  Shiny.outputBindings.register(networkOutputBinding, 'trestletech.networkbinding');
  
  </script>