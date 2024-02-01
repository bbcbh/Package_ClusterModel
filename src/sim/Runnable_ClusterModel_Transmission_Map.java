package sim;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Set;

import relationship.ContactMap;
import relationship.TransmissionMap;

public class Runnable_ClusterModel_Transmission_Map extends Runnable_ClusterModel_Transmission {

	public static final String SIM_OUTPUT_TRANMISSION_MAP = "SIM_OUTPUT_TRANMISSION_MAP";
	public static final String SIM_OUTPUT_CLUSTERS = "SIM_OUTPUT_CLUSTERS";

	protected TransmissionMap transmissionMap = null;

	public Runnable_ClusterModel_Transmission_Map(long cMap_seed, long sim_seed, int[] POP_COMPOSITION,
			ContactMap BASE_CONTACT_MAP, int NUM_TIME_STEPS_PER_SNAP, int NUM_SNAP) {
		super(cMap_seed, sim_seed, POP_COMPOSITION, BASE_CONTACT_MAP, NUM_TIME_STEPS_PER_SNAP, NUM_SNAP);

	}

	public void setTransmissionMap(TransmissionMap transmissionMap) {
		this.transmissionMap = transmissionMap;
	}

	public TransmissionMap getTransmissionMap() {
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

			Integer[] transEdge = new Integer[Abstract_Runnable_ClusterModel.LENGTH_TRANS_MAP_EDGE];
			transEdge[Abstract_Runnable_ClusterModel.TRANS_MAP_EDGE_INFECTIOUS] = infectious;
			transEdge[Abstract_Runnable_ClusterModel.TRANS_MAP_EDGE_SUSCEPTIBLE] = partner;
			transEdge[Abstract_Runnable_ClusterModel.TRANS_MAP_EDGE_START_TIME] = currentTime;
			transEdge[Abstract_Runnable_ClusterModel.TRANS_MAP_EDGE_ACT_INVOLVED] = actType;
			transmissionMap.addEdge(infectious, partner, transEdge);

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

			Set<TransmissionMap> clustersSet = new java.util.HashSet<>();
			clustersSet.add(transmissionMap);
			sim_output.put(SIM_OUTPUT_CLUSTERS, clustersSet);

			// Display clusters as CSV
			TransmissionMap[] clusters = clustersSet.toArray(new TransmissionMap[clustersSet.size()]);
			Arrays.sort(clusters, new Comparator<TransmissionMap>() {
				@Override
				public int compare(TransmissionMap o1, TransmissionMap o2) {
					return Integer.compare(o1.vertexSet().size(), o2.vertexSet().size());
				}
			});

			File clusterExport = baseDir; // new File(baseDir, String.format(DIRNAME_FORMAT_TRANSMISSION_CMAP,
											// this.sim_seed));
			clusterExport.mkdirs();

			try {
				
				File printFile;
				PrintWriter expWri;

				for (int cI = 0; cI < clusters.length; cI++) {
					TransmissionMap c = clusters[cI];

					printFile = new File(clusterExport, this.getRunnableId() == null? "": this.getRunnableId() +
							String.format(Simulation_ClusterModelTransmission.FILENAME_ALL_TRANSMISSION_CMAP,
									this.cMAP_SEED, this.sIM_SEED));

					expWri = new PrintWriter(printFile);
					expWri.println(c.toFullString());
					expWri.close();

				}
			} catch (IOException ex) {
				ex.printStackTrace(System.err);
				System.out.println("Index case:");
				System.out.println(seedInfectedStr.toString());

				for (int cI = 0; cI < clusters.length; cI++) {
					TransmissionMap c = clusters[cI];
					System.out.println(String.format("Transmission map <%d, %d>", this.sIM_SEED, cI));
					System.out.println(c.toFullString());
					System.out.println();
				}

			}
		}
	}

}
