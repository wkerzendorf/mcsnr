from mpld3 import plugins, utils
import matplotlib

css = """
table
{
  border-collapse: collapse;
}
th
{
  color: #ffffff;
  background-color: #000000;
}
td
{
  background-color: #cccccc;
}
table, th, td
{
  font-family:Arial, Helvetica, sans-serif;
  border: 1px solid black;
  text-align: right;
}
"""

class LinkedView(plugins.PluginBase):
    """A simple plugin showing how multiple axes can be linked"""

    JAVASCRIPT = """
    var LinkedViewPlugin = function(fig, prop){
      this.fig = fig;
      this.prop = mpld3.process_props(this, prop, {},
                                      ["idpts", "idline", "data"]);
    }

    LinkedViewPlugin.prototype.draw = function(){
      var pts = mpld3.get_element(this.prop.idpts);
      var line = mpld3.get_element(this.prop.idline);
      var data = this.prop.data;


      function mouseover(d, i){
        line.data = data[i];
        line.ax.set_axlim(d3.extent(data[i], function(d){return d[0];}),
                          d3.extent(data[i], function(d){return d[1];}));
        line.elements().transition()
            .attr("d", line.datafunc(line.data))
            .style("stroke", this.style.fill);
      }

      pts.elements().on("click", mouseover);
    };

    mpld3.register_plugin("linkedview", LinkedViewPlugin);
    """

    def __init__(self, points, line, linedata):
        if isinstance(points, matplotlib.lines.Line2D):
            suffix = "pts"
        else:
            suffix = None

        self.dict_ = {"type": "linkedview",
                      "idpts": utils.get_id(points, suffix),
                      "idline": utils.get_id(line),
                      "data": linedata}

