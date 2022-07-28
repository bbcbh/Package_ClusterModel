package sim;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Set;

import relationship.ContactMap;

public class Runnable_ClusterModel_Transmission_ContactMap extends Runnable_ClusterModel_Transmission {

	public static final String FILENAME_FORMAT_TRANSMISSION_CMAP = "Seed_%s_TransmissionMap_%d.csv";
	public static final String DIRNAME_FORMAT_TRANSMISSION_CMAP = "TransMap_%d";

	public static final int TRANSMAP_EDGE_INFECTIOUS = 0;
	public static final int TRANSMAP_EDGE_SUSCEPTIBLE = TRANSMAP_EDGE_INFECTIOUS + 1;
	public static final int TRANSMAP_EDGE_START_TIME = TRANSMAP_EDGE_SUSCEPTIBLE + 1;
	public static final int TRANSMAP_EDGE_ACT_INVOLVED = TRANSMAP_EDGE_START_TIME + 1;
	public static final int LENGTH_TRANSMAP_EDGE = TRANSMAP_EDGE_ACT_INVOLVED + 1;

	public static final String SIM_OUTPUT_TRANMISSION_MAP = "SIM_OUTPUT_TRANMISSION_MAP";
	public static final String SIM_OUTPUT_CLUSTERS = "SIM_OUTPUT_CLUSTERS";

	protected ContactMap transmissionMap = null;

	public Runnable_ClusterModel_Transmission_ContactMap(long seed, int[] POP_COMPOSITION, ContactMap BASE_CONTACT_MAP,
			int NUM_TIME_STEPS_PER_SNAP, int SNAP_FREQ) {
		super(seed, POP_COMPOSITION, BASE_CONTACT_MAP, NUM_TIME_STEPS_PER_SNAP, SNAP_FREQ);

	}

	public void setTransmissionMap(ContactMap transmissionMap) {
		this.transmissionMap = transmissionMap;
	}

	public ContactMap getTransmissionMap() {
		return transmissionMap;
	}

	@Override
	protected void transmission_success(int currentTime, Integer infectious, int partner, int site_target, int actType,
			Object[] simulation_store) {

		super.transmission_success(currentTime, infectious, partner, site_target, actType, simulation_store);

		if (transmissionMap != null) {

			if (!transmissionMap.containsVertex(infectious)) {
				transmissionMap.addVertex(infectious);
			}
			if (!transmissionMap.containsVertex(partner)) {
				transmissionMap.addVertex(partner);
			}

			Integer[] existEdge = transmissionMap.getEdge(infectious, partner);

			if (existEdge == null) {
				existEdge = new Integer[LENGTH_TRANSMAP_EDGE];
				existEdge[TRANSMAP_EDGE_INFECTIOUS] = infectious;
				existEdge[TRANSMAP_EDGE_SUSCEPTIBLE] = partner;
				existEdge[TRANSMAP_EDGE_START_TIME] = currentTime;
				existEdge[TRANSMAP_EDGE_ACT_INVOLVED] = 1 << actType;
				transmissionMap.addEdge(infectious, partner, existEdge);
			} else {
				existEdge[TRANSMAP_EDGE_ACT_INVOLVED] |= 1 << actType;
			}
		}
	}

	@Override
	protected Object[] preSimulation() {
		Object[] simulation_store = super.preSimulation();
		int offset = simulation_store.length;
		simulation_store = Arrays.copyOf(simulation_store, simulation_store.length + 1);
		int[][] seedInfected = new int[currently_infectious.length][];
		// Store initially infected
		for (int site = 0; site < currently_infectious.length; site++) {
			ArrayList<Integer> currently_infectious_by_site = currently_infectious[site];
			seedInfected[site] = new int[currently_infectious_by_site.size()];
			int c = 0;
			for (Integer infectious : currently_infectious_by_site) {
				seedInfected[site][c] = infectious;
				c++;
			}
		}
		simulation_store[offset] = seedInfected;
		return simulation_store;
	}

	@Override
	protected void postSimulation(Object[] simulations_store) {
		super.postSimulation(simulations_store);

		// Store index case(s)
		StringBuilder seedInfectedStr = new StringBuilder();
		int[][] seedInfected = (int[][]) simulations_store[0];
		for (int site = 0; site < seedInfected.length; site++) {
			for (int i = 0; i < seedInfected[site].length; i++) {
				seedInfectedStr.append(site);
				seedInfectedStr.append(',');
				seedInfectedStr.append(seedInfected[site][i]);
				seedInfectedStr.append('\n');
			}
		}

		if (transmissionMap != null) {

			sim_output.put(SIM_OUTPUT_TRANMISSION_MAP, transmissionMap);
			Set<ContactMap> clustersSet = transmissionMap.getContactCluster();
			sim_output.put(SIM_OUTPUT_CLUSTERS, clustersSet);

			// Display clusters as CSV
			ContactMap[] clusters = clustersSet.toArray(new ContactMap[clustersSet.size()]);
			Arrays.sort(clusters, new Comparator<ContactMap>() {
				@Override
				public int compare(ContactMap o1, ContactMap o2) {
					return Integer.compare(o1.vertexSet().size(), o2.vertexSet().size());
				}
			});

			File clusterExport = new File(baseDir, String.format(DIRNAME_FORMAT_TRANSMISSION_CMAP, this.seed));
			clusterExport.mkdirs();

			File printFile;
			PrintWriter expWri;

			try {

				printFile = new File(clusterExport, String.format(FILENAME_FORMAT_INDEX_CASE_LIST, this.seed));

				expWri = new PrintWriter(printFile);
				expWri.println(seedInfectedStr.toString());
				expWri.close();

				for (int cI = 0; cI < clusters.length; cI++) {
					ContactMap c = clusters[cI];

					printFile = new File(clusterExport,
							String.format(FILENAME_FORMAT_TRANSMISSION_CMAP, Long.toString(this.seed), cI));

					expWri = new PrintWriter(printFile);
					expWri.println(c.toFullString());
					expWri.close();

				}
			} catch (IOException ex) {
				ex.printStackTrace(System.err);
				System.out.println("Index case:");
				System.out.println(seedInfectedStr.toString());

				for (int cI = 0; cI < clusters.length; cI++) {
					ContactMap c = clusters[cI];
					System.out.println(String.format("Transmission map <%d, %d>", this.seed, cI));
					System.out.println(c.toFullString());
					System.out.println();
				}

			}
		}
	}

}
