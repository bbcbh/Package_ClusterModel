package visualise;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import javax.imageio.ImageIO;
import javax.swing.JFrame;

import org.jgrapht.ext.JGraphXAdapter;

import com.mxgraph.layout.mxFastOrganicLayout;
import com.mxgraph.layout.mxIGraphLayout;
import com.mxgraph.model.mxICell;
import com.mxgraph.swing.mxGraphComponent;
import com.mxgraph.util.mxCellRenderer;
import com.mxgraph.util.mxConstants;
import com.mxgraph.util.mxStyleUtils;

import relationship.ContactMap;

public class Visualise_ContactMap extends JGraphXAdapter<Integer, Integer[]> {

	protected static final String EDGE_STYLE = "EDGE_STYLE";

	public Visualise_ContactMap(ContactMap cmap) {
		super(cmap);
		Map<String, Object> edge_style = this.getStylesheet().getDefaultEdgeStyle();
		edge_style.put(mxConstants.STYLE_FONTSIZE, 0);
		// edge_style.put(mxConstants.STYLE_ENDARROW, mxConstants.ARROW_OVAL);
		// edge_style.put(mxConstants.STYLE_STARTARROW, mxConstants.ARROW_OVAL);
		edge_style.put(mxConstants.STYLE_STARTSIZE, 0);
		edge_style.put(mxConstants.STYLE_ENDSIZE, 0);

	}
	
	public static final int FLAG_GEN_CSV = 1 << 0;
	public static final int FLAG_DISP_FIG = 1 << 1;
	public static final int FLAG_GEN_PNG = 1 << 2;

	public static void visualiseCluster(File baseDir, int[] popComposition, ContactMap cluster, Long seed,
			int clusterCount, int flag) throws IOException {
		
		if ((flag & FLAG_GEN_CSV) > 0) {

			String mapRep = cluster.toFullString();
			BufferedWriter fileWri = new BufferedWriter(new FileWriter(new File(baseDir,
					String.format("cluster_%d_%d_%d.csv", seed, clusterCount, cluster.vertexSet().size()))));
			fileWri.append(mapRep);
			fileWri.close();
		}

		if (((flag & FLAG_DISP_FIG) > 0) || ((flag & FLAG_GEN_PNG) > 0)) {

			Visualise_ContactMap mxGraph = new Visualise_ContactMap(cluster);
			HashMap<Integer, mxICell> vMap = mxGraph.getVertexToCellMap();

			for (Integer v : vMap.keySet()) {
				mxICell cell = vMap.get(v);

				if (v >= popComposition[0] + popComposition[1] + popComposition[2]) {
					mxStyleUtils.setCellStyles(mxGraph.getModel(), new mxICell[] { cell }, mxConstants.STYLE_FILLCOLOR,
							"#00FFFF");					
				} else if (v >= popComposition[0] + popComposition[1]) {
					mxStyleUtils.setCellStyles(mxGraph.getModel(), new mxICell[] { cell }, mxConstants.STYLE_FILLCOLOR,
							"#00FF00");
				} else if (v < popComposition[0]) {
					mxStyleUtils.setCellStyles(mxGraph.getModel(), new mxICell[] { cell }, mxConstants.STYLE_FILLCOLOR,
							"#FF0000");
				}

			}

			mxIGraphLayout layout = new mxFastOrganicLayout(mxGraph);
			layout.execute(mxGraph.getDefaultParent());
			mxGraph.setCellsLocked(true);

			if ((flag & FLAG_DISP_FIG) > 0) {

				if (cluster.vertexSet().size() < 1000) {

					mxGraphComponent graphComponent = new mxGraphComponent(mxGraph);
					JFrame dispFrame = new JFrame(
							String.format("Cluser size = %d Seed #%d: ", cluster.vertexSet().size(), seed));
					dispFrame.getContentPane().add(graphComponent);

					dispFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
					dispFrame.setSize(400, 320);
					dispFrame.setVisible(true);
				}
			}

			if ((flag & FLAG_GEN_PNG) > 0) {
				BufferedImage img = mxCellRenderer.createBufferedImage(mxGraph, null, 1, Color.WHITE, false, null);
				ImageIO.write(img, "PNG", new File(baseDir, String.format("Cluster_%d_%d.png", seed, clusterCount)));
			}
		}
	}

	/*
	 * public static void main(String[] arg) { ContactMap test = new ContactMap();
	 * for (int i = 0; i < 10; i++) { test.addVertex(i + 1); }
	 * 
	 * test.addEdge(1, 2); test.addEdge(2, 3); test.addEdge(3, 1); test.addEdge(4,
	 * 5); test.addEdge(6, 8); test.addEdge(8, 9);
	 * 
	 * Visualise_ContactMap mxGraph = new Visualise_ContactMap(test);
	 * 
	 * HashMap<Integer, mxICell> vMap = mxGraph.getVertexToCellMap(); for (Integer v
	 * : vMap.keySet()) { mxICell cell = vMap.get(v); if ((v % 2) == 0) {
	 * 
	 * mxStyleUtils.setCellStyles(mxGraph.getModel(),new mxICell[] {cell},
	 * mxConstants.STYLE_FILLCOLOR, "#FF0000"); } }
	 * 
	 * mxIGraphLayout layout = new mxCircleLayout(mxGraph);
	 * layout.execute(mxGraph.getDefaultParent()); mxGraph.setCellsLocked(true);
	 * 
	 * 
	 * mxGraphComponent graphComponent = new mxGraphComponent(mxGraph); JFrame
	 * testFrame = new JFrame("Test");
	 * testFrame.getContentPane().add(graphComponent);
	 * 
	 * testFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	 * testFrame.setSize(400, 320); testFrame.setVisible(true); }
	 */
}
