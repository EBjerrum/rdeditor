// make sure our module is only defined
// only once.
require.undef('molsvg_widget');

// Define the `molsvg_widget` module using the Jupyter widgets framework.
define('molsvg_widget', ["@jupyter-widgets/base"],
       function(widgets) {

    // The frontend class:
    var MolSVGView = widgets.DOMWidgetView.extend({

        // This method creates the HTML widget.
        render: function() {
            this.svg_div = document.createElement('div');
            this.el.appendChild(this.svg_div);
            this.model.on('change:svg', this.svg_changed, this);
            this.svg_changed();
        },
        
        // called when the SVG is updated on the Python side
        svg_changed: function() {
            var txt = this.model.get('svg'); 
            this.svg_div.innerHTML = txt;
            var sels = this.svg_div.getElementsByClassName("atom-selector");
            for(var i=0;i<sels.length;i++){
                sels[i].onclick = (evt) => { return this.atom_clicked(evt) };
            }
            var sels = this.svg_div.getElementsByClassName("bond-selector");
            for(var i=0;i<sels.length;i++){
                sels[i].onclick = (evt) => { return this.bond_clicked(evt) };
            }
            
        },

        // callback for when an atom is clicked
        atom_clicked: function(evt) {
            //alert("  "+evt+"|"+this);
            if(!evt.currentTarget.getAttribute('class')){
                return;
            }
            var satmid = evt.currentTarget.getAttribute('class').match(/atom-([0-9]+)/);
            if(satmid.length >1){
                var atmid = Number(satmid[1]);
                this.model.set('clicked_atom_idx','Event');
                this.touch();
                this.model.set('clicked_atom_idx',String(atmid));
                this.touch();
            }
         },
        // callback for when a bond is clicked        
        bond_clicked: function(evt) {
            //alert("  "+evt+"|"+this);
            if(!evt.currentTarget.getAttribute('class')){
                return;
            }
            var sbondid = evt.currentTarget.getAttribute('class').match(/bond-([0-9]+)/);
            if(sbondid.length >1){
                var bondid = Number(sbondid[1]);
                this.model.set('clicked_bond_idx','Event');
                this.touch();
                this.model.set('clicked_bond_idx',String(bondid));
                this.touch();
            }
        },


    });

    return {
        MolSVGView : MolSVGView
    };
});