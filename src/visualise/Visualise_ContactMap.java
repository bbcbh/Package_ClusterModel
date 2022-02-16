package visualise;

import java.util.Map;

import org.jgrapht.ext.JGraphXAdapter;

import com.mxgraph.util.mxConstants;

import relationship.ContactMap;

public class Visualise_ContactMap extends JGraphXAdapter<Integer, Integer[]> {

	protected static final String EDGE_STYLE = "EDGE_STYLE";

	public Visualise_ContactMap(ContactMap cmap) {
		super(cmap);

		Map<String, Object> edge_style = this.getStylesheet().getDefaultEdgeStyle();
		edge_style.put(mxConstants.STYLE_FONTSIZE, 0);
		//edge_style.put(mxConstants.STYLE_ENDARROW, mxConstants.ARROW_OVAL);
		//edge_style.put(mxConstants.STYLE_STARTARROW, mxConstants.ARROW_OVAL);
		edge_style.put(mxConstants.STYLE_STARTSIZE, 0);
		edge_style.put(mxConstants.STYLE_ENDSIZE, 0);

	}
/*
	public static void main(String[] arg) {
		ContactMap test = new ContactMap();
		for (int i = 0; i < 10; i++) {
			test.addVertex(i + 1);
		}

		test.addEdge(1, 2);
		test.addEdge(2, 3);
		test.addEdge(3, 1);
		test.addEdge(4, 5);
		test.addEdge(6, 8);
		test.addEdge(8, 9);

		Visualise_ContactMap mxGraph = new Visualise_ContactMap(test);		
			
		HashMap<Integer, mxICell> vMap = mxGraph.getVertexToCellMap();
		for (Integer v : vMap.keySet()) {
			mxICell cell = vMap.get(v);		
			if ((v % 2) == 0) {									
				
				mxStyleUtils.setCellStyles(mxGraph.getModel(),new mxICell[] {cell}, mxConstants.STYLE_FILLCOLOR, "#FF0000");
			}
		}		

		mxIGraphLayout layout = new mxCircleLayout(mxGraph);
		layout.execute(mxGraph.getDefaultParent());		
		mxGraph.setCellsLocked(true);
		
		
		mxGraphComponent graphComponent = new mxGraphComponent(mxGraph);
		JFrame testFrame = new JFrame("Test");
		testFrame.getContentPane().add(graphComponent);

		testFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		testFrame.setSize(400, 320);
		testFrame.setVisible(true);
	}
*/
}
