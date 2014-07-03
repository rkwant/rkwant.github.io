// JavaScript Document
$(document).ready(function () {

    
    $("#addbutton").click(function () {
        var counter = $("#species tr").length;
        var newrow = $("<tr>");

        var cols = "<td><input type='text' value='1' class='mult' size='4'></input></td>";
        cols += "<td><input type='text' value='3.5' class='ka' size='4'></input></td>";
        cols += "<td><input type='text' value='3.5' class='ks' size='4'></input></td>";
        cols += "<td><input type='text' value='0.00001' class='conc' size='4'></input></td>";
        cols += "<td><a href='#' class='remover'>Delete</a></td>";
        newrow.append(cols);
        $("#species").append(newrow);
    });

    $("#addgradient").click(function () {
        var counter = $("#gradient tr").length;
        var newrow = $("<tr>");

        var cols = "<td><input type='text' value='1' class='start' size='4'></input></td>";
        cols += "<td><input type='text' value='1' class='duration' size='4'></input></td>";
        cols += "<td><input type='text' value='0.01' class='endconc' size='4'></input></td>";
        cols += "<td><a href='#' class='remover'>Delete</a></td>";
        newrow.append(cols);
        $("#gradient").append(newrow);
    });

    $("#gradient").on("click", "tr a", function (e) {
        $(this).parents("tr").remove();
    });

    $("#species").on("click", "tr a", function (e) {
        $(this).parents("tr").remove();
    });
	
	
	
	function get_gradient_addition(cv, ts, pf, gsl) {
        var amount = 0;
		var av,bv,concv;
		
        for (var k=0; k<gsl.length; k++)
		{
            av = gsl[k][0];
            bv = gsl[k][1];//va[1];
            concv=gsl[k][2]; //va[2];
            amount = amount + pf*ts*concv*((0.5*(cv-av)+0.5*Math.abs(cv-av))-(0.5*(cv-(av+bv))+0.5*Math.abs(cv-(av+bv))))/(bv)
        }
    //$("#tst").html(""+av+" "+bv+" "+concv);
        return amount;
    }
    

    
    

    function get_sol_conc(f, Fx, kax, Fa, kas, ms) {
        var term1 = numeric.mul(Math.pow(2, -f), Fa);
        var term2 = numeric.div(1, numeric.mul(f, kax, Fa));
        var term3 = numeric.add(-1, numeric.mul(f, kax, Fa), -kax * Fx, numeric.mul(-kas, ms));
        var term4 = numeric.sqrt(numeric.add(numeric.mul(-4, f, kax, Fa, numeric.add(-1, numeric.mul(-kas, ms))), numeric.pow(numeric.add(1, numeric.mul(-1, f, kax, Fa), kax * Fx, numeric.mul(kas, ms)), 2)));
        var result = numeric.mul(term1, numeric.pow(numeric.mul(term2, numeric.add(term3, term4)), f));

        for (var j = 0; j < result.length; j++) {
            if (isNaN(result[j])) result[j] = 0;
        }

        return result;
    }

    function run_simulation(elution_vol, col_vol, multiplicity, plates, resin, kax, kas, protein_start, timestep, gsl) {
        var vol_per_plate = col_vol / plates;
        var vs = elution_vol;
        var flow = 0.2;
        var plate_flow = flow / vol_per_plate;

        var iterations = Math.floor(vs / (timestep * flow));

        var volumes = new Array(45);
        for (var i = 0; i < iterations; i++) {
            volumes[i] = i;
        }
        volumes = numeric.mul(volumes, flow, timestep);

        var solution = new Array(plates);
        var bound = new Array(plates);
        var eluent = new Array(plates);
        for (i = 0; i < plates; i++) {
            solution[i] = 0;
            bound[i] = 0;
            eluent[i] = 0;
        }

        solution[0] = protein_start;

        var out = new Array();
		var gradout = new Array();
		
		var tr = {
			chrom: {
            	data: []
        	}, 
			grad:{
				data:[]
			}
		};
		
        var current_volume = 0;
        var total = [];
        var solution_change = [];
        var eluent_change = [];

        for (i = 0; i < iterations; i++) {
            current_volume = i * timestep * flow;
            total = numeric.add(solution, bound);
            solution = get_sol_conc(multiplicity, resin, kax, total, kas, eluent);
            bound = numeric.sub(total, solution)

            solution_change = numeric.mul(plate_flow, solution, timestep);
            eluent_change = numeric.mul(plate_flow, eluent, timestep);

            solution = numeric.sub(solution, solution_change);
            tr.chrom.data.push({
                x: current_volume,
                y: solution_change.pop()
            });
            solution_change.unshift(0);
            solution = numeric.add(solution, solution_change);

            eluent = numeric.sub(eluent, eluent_change);
            tr.grad.data.push({
                x: current_volume,
                y: eluent_change.pop() / (plate_flow*timestep)
            });
            eluent_change.unshift(0);
            eluent = numeric.add(eluent, eluent_change);
            eluent[0] = eluent[0] + get_gradient_addition(current_volume, timestep, plate_flow, gsl);
            
        }
        return tr;
    }






    $("#plot").click(function () {
        //define the data object to be plotted
		while (chromatoData.length > 0) chromatoData.pop();
		gradData.pop();
		
        //get the gradient information
        var gradient_segment_list = []
        $("#gradient tr").each(function (ind, value) {
            if (ind > 0) {
                gradient_segment_list.push([
                parseFloat($(value).find(".start").val()),
                parseFloat($(value).find(".duration").val()),
                parseFloat($(value).find(".endconc").val())])
            }
            //console.log(gradient_segment_list);
        });
		
		var d = new Date();
		var n = d.getTime();
        $("#species tr").each(function (ind, value) {

            if (ind > 0) {
                //console.log($(value).find(".mult").val());

                var cur_data = run_simulation(
                parseFloat($("#elution_volume").val()),
				parseFloat($("#volume").val()),
                parseFloat($(value).find(".mult").val()),
                parseFloat($("#plates").val()),
                parseFloat($("#resin").val()),
                Math.pow(10, parseFloat($(value).find(".ka").val())),
                Math.pow(10, parseFloat($(value).find(".ks").val())),
                parseFloat($(value).find(".conc").val()),
                parseFloat($("#timestep").val()),
				gradient_segment_list);
				
				chromatoData.push({data: cur_data.chrom.data, color: "#c05020"});
				if (ind==1) gradData.push({data: cur_data.grad.data, color: "#c05020"});
            }
        });
		var d = new Date();
		$("#tst").html("Completed in "+(d.getTime() - n)+" ms");
		chromatograph.update();
		gradgraph.update();
        //$("#chartContainer").dxChart(toPlot);
		//$("#gradContainer").dxChart(gradPlot);
		
    });
    
	var chromatoData = [
        {
            data: [ { x: 0, y: 120 }, { x: 1, y: 890 }, { x: 2, y: 38 }, { x: 3, y: 70 }, { x: 4, y: 32 } ],
            color: "#c05020"
        }
	];
	
	var gradData = [
        {
            data: [ { x: 0, y: 120 }, { x: 1, y: 890 }, { x: 2, y: 38 }, { x: 3, y: 70 }, { x: 4, y: 32 } ],
            color: "#c05020"
        }
	];


// instantiate our graph!

var chromatograph = new Rickshaw.Graph( {
	element: document.getElementById("chartContainer"),
	renderer: 'line',
	series: chromatoData
} );

var gradgraph = new Rickshaw.Graph( {
	element: document.getElementById("gradContainer"),
	renderer: 'line',
	series: gradData
} );

var x_ticks = new Rickshaw.Graph.Axis.X( {
	graph: chromatograph,
	orientation: 'bottom',
	element: document.getElementById('xaxis'),
	pixelsPerTick: 200,
	tickFormat: Rickshaw.Fixtures.Number.formatKMBT
} );



chromatograph.render();
gradgraph.render();


    $("#addbutton").click();
	$("#plot").click();

});