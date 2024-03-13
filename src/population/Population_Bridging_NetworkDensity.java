package population;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;

import person.AbstractIndividualInterface;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;
import sim.Abstract_Runnable_ClusterModel;

public class Population_Bridging_NetworkDensity extends Population_Bridging_Scheduled {

	private static final long serialVersionUID = -3797890258608814830L;

	// As defined in Karang et al. doi:10.53638/phpma.2017.v5.i1.p15
	private static final int NETWORK_DENSITY_LVL_3 = 3; // at least 3 sexual relations,
	private static final int NETWORK_DENSITY_LVL_4 = NETWORK_DENSITY_LVL_3 + 1; // at least 4 sexual relations with one
																				// common partner
	private static final int NETWORK_DENSITY_LVL_5 = NETWORK_DENSITY_LVL_4 + 1; // at least 5 sexual relations and 2
																				// common partners
	private static final int NETWORK_DENSITY_LVL_6 = NETWORK_DENSITY_LVL_5 + 1; // at least 6 sexual relations where at
																				// least 2 common partners
	// also has partnership

	private static final int[] N_PER_LVL = new int[] { 99, 20, 5, 1 };
	private static float N_SUM;

	private RandomGenerator RNG_density;

	ContactMap cMapLast12Months = new ContactMap();

	public Population_Bridging_NetworkDensity(long seed) {
		super(seed);

		float nSum = 0;
		for (int n : N_PER_LVL) {
			nSum += n;
		}
		N_SUM = nSum;
		RNG_density = new MersenneTwisterRandomGenerator(this.getSeed());

	}

	@Override
	protected void checkContactMaps(Integer[] link, ContactMap[] cMaps) {
		super.checkContactMaps(link, cMaps);

		if (link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] > 0) {

			int link_cMap_expiry_time = link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]
					+ link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION]
					+ AbstractIndividualInterface.ONE_YEAR_INT;
			for (int i = 0; i < 2; i++) {
				if (!cMapLast12Months.containsVertex(link[i])) {
					cMapLast12Months.addVertex(link[i]);
				}
			}
			boolean contain_edge = cMapLast12Months.containsEdge(link[0], link[1]);

			if (contain_edge) {
				Set<Integer[]> existingEdge = cMapLast12Months.getAllEdges(link[0], link[1]);
				boolean firstEdge = true;
				for (Integer[] e : existingEdge) {
					if (firstEdge) {
						// Replace same edge with longer expiry time
						int org_cMap_expiry_time = e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]
								+ e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION]
								+ AbstractIndividualInterface.ONE_YEAR_INT;
						if (link_cMap_expiry_time > org_cMap_expiry_time) {
							e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME] = link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME];
							e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] = link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION];
						}
						firstEdge = false;
					} else {
						// Should only be one edge anyway
						cMapLast12Months.removeEdge(e);
					}
				}
			} else {
				cMapLast12Months.addEdge(link[0], link[1], link);
			}
		}

	}

	@Override
	public void advanceTimeStep(int deltaT) {
		super.advanceTimeStep(deltaT);

		boolean reportPartnerStat = lastPartnershipScheduling < AbstractIndividualInterface.ONE_YEAR_INT
				? (getGlobalTime() == AbstractIndividualInterface.ONE_YEAR_INT)
				: getGlobalTime() == lastPartnershipScheduling;

		// All Key = PID
		HashMap<Integer, Integer[]> lookup_partners = new HashMap<>();
		HashMap<Integer, Integer> lookup_densityLevel = new HashMap<>();
		HashMap<Integer, Integer> lookup_firstMatchFrom = new HashMap<>();
		HashMap<Integer, ArrayList<Integer[]>> lookup_no_partnership_between_partners = new HashMap<>();

		// Clear out all partner that finished 12 months ago
		ArrayList<Integer[]> edge_to_remove = new ArrayList<>();

		for (Integer[] e : cMapLast12Months.edgeSet()) {
			if (getGlobalTime() > AbstractIndividualInterface.ONE_YEAR_INT
					&& e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]
							+ e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] < getGlobalTime()
									- AbstractIndividualInterface.ONE_YEAR_INT) {
				edge_to_remove.add(e);
			}
		}

		for (Integer[] e : edge_to_remove) {
			cMapLast12Months.removeEdge(e);
		}

		for (Integer index_pid : cMapLast12Months.vertexSet()) {

			if (cMapLast12Months.degreeOf(index_pid) >= 3) {

				ArrayList<Integer> connection = new ArrayList<>();
				Set<Integer[]> edgeAll = cMapLast12Months.edgesOf(index_pid);

				int firstMatchAt = Integer.MAX_VALUE;
				int validUntil = Integer.MAX_VALUE;
				int density_lvl = NETWORK_DENSITY_LVL_3;

				for (Integer[] e : edgeAll) {
					Integer partner = e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1].equals(index_pid)
							? e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2]
							: e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1];

					int pt = Collections.binarySearch(connection, partner);
					if (pt < 0) {
						connection.add(~pt, partner);
					}
					firstMatchAt = Math.min(firstMatchAt,
							e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]);
					validUntil = Math.min(validUntil,
							e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]
									+ e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION]
									+ AbstractIndividualInterface.ONE_YEAR_INT);
				}

				Integer[] partners = connection.toArray(new Integer[0]);
				ArrayList<Integer[]> no_partnership_array_index = new ArrayList<>();

				for (int p1 = 0; p1 < partners.length; p1++) {
					for (int p2 = p1 + 1; p2 < partners.length; p2++) {
						if (!cMapLast12Months.containsEdge(partners[p1], partners[p2])) {
							no_partnership_array_index.add(new Integer[] { partners[p1], partners[p2] });
						} else {
							density_lvl = Math.min(density_lvl + 1, NETWORK_DENSITY_LVL_6);
						}
					}
				}
				lookup_partners.put(index_pid, partners);
				lookup_no_partnership_between_partners.put(index_pid, no_partnership_array_index);
				lookup_densityLevel.put(index_pid, density_lvl);
				lookup_firstMatchFrom.put(index_pid, firstMatchAt);

			}

		}

		// Analyse density from individual level
		if (lookup_densityLevel.size() > N_SUM) {

			// Calculate current and ideal density count
			int[] count = new int[N_PER_LVL.length];
			for (Integer lvl : lookup_densityLevel.values()) {
				count[lvl.intValue() - NETWORK_DENSITY_LVL_3]++;
			}
			int[] idealCount = new int[N_PER_LVL.length];
			for (int i = 0; i < idealCount.length; i++) {
				idealCount[i] = Math.round(lookup_densityLevel.size() * N_PER_LVL[i] / N_SUM);
			}

			// Allocated candidate, order by number of connections
			ArrayList<ArrayList<Integer>> densityLvlCandidates = new ArrayList<>(N_PER_LVL.length);
			for (int i = 0; i < N_PER_LVL.length; i++) {
				densityLvlCandidates.add(new ArrayList<>());
			}
			for (Integer index_pid : lookup_densityLevel.keySet()) {
				int lvl_pt = lookup_densityLevel.get(index_pid).intValue() - NETWORK_DENSITY_LVL_3;
				ArrayList<Integer> index_arr = densityLvlCandidates.get(lvl_pt);
				index_arr.add(index_pid);
			}

			String init_count = Arrays.toString(count);
			for (int dI = N_PER_LVL.length - 1; dI > 0; dI--) {

				while (idealCount[dI] > count[dI]) {
					int colSel = 0;
					while (colSel < densityLvlCandidates.size() && densityLvlCandidates.get(colSel).size() == 0) {
						colSel++;
					}

					if (colSel >= densityLvlCandidates.size()) {
						// No more suitable candidate
						break;
					}

					int indexSel = densityLvlCandidates.get(colSel)
							.remove(RNG_density.nextInt(densityLvlCandidates.get(colSel).size()));

					int partnerBeforeTime = lookup_firstMatchFrom.get(indexSel)
							+ AbstractIndividualInterface.ONE_YEAR_INT;

					if (partnerBeforeTime > getGlobalTime()) {

						ArrayList<Integer[]> no_partnership_index = lookup_no_partnership_between_partners
								.get(indexSel);
						int numConnectionToAdd = (int) Math.min(dI - colSel, no_partnership_index.size());

						// Schedule a future partnership to schedule
						while (numConnectionToAdd > 0) {

							Integer[] selectPartner = no_partnership_index
									.remove(RNG_density.nextInt(no_partnership_index.size()));

							int partnershipFormTime = getGlobalTime()
									+ RNG_density.nextInt(partnerBeforeTime - getGlobalTime());
							Integer[] newEdge = new Integer[] { selectPartner[0], selectPartner[1],
									CANDIDATE_ARRAY_SOUGHT_ANY };

							ArrayList<Integer[]> schedule_partnership_ent = schedule_partnership
									.get(partnershipFormTime);

							if (schedule_partnership_ent == null) {
								schedule_partnership_ent = new ArrayList<>();
								schedule_partnership.put(partnershipFormTime, schedule_partnership_ent);
							}
							int entPt = Collections.binarySearch(schedule_partnership_ent, newEdge,
									schedule_partnership_comparator);
							if (entPt < 0) {
								schedule_partnership_ent.add(~entPt, newEdge);
							}
							numConnectionToAdd--;
						}
					}
					count[colSel]--;
					count[dI]++;
				}
			}

			if (reportPartnerStat) {
				System.out.printf("Time = %d. Unadj. density=%s. Ideal. density = %s.\n", getGlobalTime(), init_count,
						Arrays.toString(idealCount));
				if (printStatus != null) {
					if (reportPartnerStat) {
						String str_networkDensity = String.format("Network densitiy at Time %d = %s", getGlobalTime(),
								init_count);
						for (PrintStream pWri : printStatus) {
							pWri.println(str_networkDensity);
							pWri.println();
						}

					}
				}
			}

		}

	}

}
