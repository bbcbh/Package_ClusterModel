package population;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Set;

import org.apache.commons.math3.util.CombinatoricsUtils;

import person.AbstractIndividualInterface;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;
import sim.Abstract_Runnable_ClusterModel;
import util.ArrayUtilsRandomGenerator;

public class Population_Bridging_NetworkDensity extends Population_Bridging_Scheduled {

	private static final long serialVersionUID = -3797890258608814830L;

	// As defined in Karang et al. doi:10.53638/phpma.2017.v5.i1.p15
	private static final int NETWORK_DENSITY_LVL_3 = 3; // at least 3 sexual relations,
	private static final int NETWORK_DENSITY_LVL_4 = 4; // at least 4 sexual relations with one common partner
	private static final int NETWORK_DENSITY_LVL_5 = 5; // at least 5 sexual relations and 2 common partners
	private static final int NETWORK_DENSITY_LVL_6 = 6; // at least 6 sexual relations where at least 2 common partners
														// also has partnership

	private static final int[] N_PER_LVL = new int[] { 99, 20, 5, 1 };
	private static float N_SUM;


	private RandomGenerator RNG_density;

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
	public void advanceTimeStep(int deltaT) {
		super.advanceTimeStep(deltaT);

		ContactMap cMapLast12Months = contactMapInLast12Months(
				((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_ALL]);

		if (cMapLast12Months != null) {

			HashMap<Integer, ContactMap> lookup_index_connections = new HashMap<>();
			HashMap<Integer, Integer> lookup_densityLevel = new HashMap<>();
			HashMap<Integer, Integer> lookup_firstMatchFrom = new HashMap<>();

			for (Integer pId : cMapLast12Months.vertexSet()) {
				if (cMapLast12Months.degreeOf(pId) >= 3) {
					Set<Integer[]> edgeAll = cMapLast12Months.edgesOf(pId);
					if (edgeAll.size() >= 3) {
						ContactMap indivdualMap = new ContactMap();
						int firstMatchAt = Integer.MAX_VALUE;
						for (Integer[] e : edgeAll) {
							firstMatchAt = Math.min(firstMatchAt,
									e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]);

							if (!indivdualMap.containsVertex(e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1])) {
								indivdualMap.addVertex(e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1]);
							}
							if (!indivdualMap.containsVertex(e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2])) {
								indivdualMap.addVertex(e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2]);
							}

							indivdualMap.addEdge(e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1],
									e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2], e);
						}
						lookup_index_connections.put(pId, indivdualMap);
						lookup_firstMatchFrom.put(pId, firstMatchAt);
					}

					Integer[] indexCaseArr = lookup_index_connections.keySet()
							.toArray(new Integer[lookup_index_connections.size()]);

					for (Integer indexCase : indexCaseArr) {
						ContactMap indivdualMap = lookup_index_connections.get(indexCase);
						Set<Integer> vertices = indivdualMap.vertexSet();
						int density_lvl = NETWORK_DENSITY_LVL_3;

						for (Integer v : vertices) {
							if (!v.equals(indexCase)) {
								ArrayList<Integer> third_partner_list = new ArrayList<>();

								// Check how many partner of partner is also a direct partner to index case
								Set<Integer[]> vEdgeAll = cMapLast12Months.edgesOf(v);
								for (Integer[] e : vEdgeAll) {
									Integer third_partner = v
											.equals(e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1])
													? e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2]
													: e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1];
									if (!third_partner.equals(indexCase)) {
										boolean isCommonPartner = hasEdgeIn12Months(cMapLast12Months, indexCase,
												third_partner);
										if (isCommonPartner) {
											third_partner_list.add(third_partner);
											density_lvl = Math.max(density_lvl, NETWORK_DENSITY_LVL_4);
										}
									}

								}
								if (third_partner_list.size() > 1) {
									density_lvl = Math.max(density_lvl, NETWORK_DENSITY_LVL_5);
									// Check each partner combination and see if they are partner as well
									for (int i = 0; i < third_partner_list.size()
											&& density_lvl < NETWORK_DENSITY_LVL_6; i++) {
										for (int j = i; j < third_partner_list.size()
												&& density_lvl < NETWORK_DENSITY_LVL_6; j++) {
											if (hasEdgeIn12Months(cMapLast12Months, third_partner_list.get(i),
													third_partner_list.get(j))) {
												density_lvl = Math.max(density_lvl, NETWORK_DENSITY_LVL_6);
											}
										}
									}
								}
							}
						}
						lookup_densityLevel.put(indexCase, density_lvl);
					}
				}
			}

			// Analyse density from individual level
			if (lookup_densityLevel.size() > N_SUM) {
				int[] count = new int[N_PER_LVL.length];
				for (Integer lvl : lookup_densityLevel.values()) {
					count[lvl.intValue() - NETWORK_DENSITY_LVL_3]++;
				}

				int[] idealCount = new int[N_PER_LVL.length];
				for (int i = 0; i < idealCount.length; i++) {
					idealCount[i] = Math.round(lookup_densityLevel.size() * N_PER_LVL[i] / N_SUM);
				}
				ArrayList<ArrayList<Integer>> densityLvlCandidates = new ArrayList<>(N_PER_LVL.length);
				for (int i = 0; i < N_PER_LVL.length; i++) {
					densityLvlCandidates.add(new ArrayList<>());
				}
				for (Integer indexId : lookup_densityLevel.keySet()) {
					int lvl_pt = lookup_densityLevel.get(indexId).intValue() - NETWORK_DENSITY_LVL_3;
					ArrayList<Integer> index_arr = densityLvlCandidates.get(lvl_pt);
					int arrPt = Collections.binarySearch(index_arr, indexId, new Comparator<Integer>() {
						@Override
						public int compare(Integer o1, Integer o2) {
							int res = -Integer.compare(cMapLast12Months.degreeOf(o1), cMapLast12Months.degreeOf(o2));
							if (res == 0) {
								res = Integer.compare(o1, o2);
							}
							return res;
						}
					});
					if (arrPt < 0) {
						index_arr.add(~arrPt, indexId);
					}
				}
				
				String init_count = Arrays.toString(count);
				for (int dI = N_PER_LVL.length - 1; dI > 0; dI--) {
					
					while (idealCount[dI] > count[dI]) {
						int colSel = 0;
						while (colSel < densityLvlCandidates.size() && densityLvlCandidates.get(colSel).size() == 0) {
							colSel++;
						}
						
						if(colSel >= densityLvlCandidates.size()) {
							// No more suitable candidate
							break;
						}
						
						int indexSel = densityLvlCandidates.get(colSel).remove(0);
						ContactMap indivdual_map = lookup_index_connections.get(indexSel);

						ArrayList<Integer> partner_list = new ArrayList<>(indivdual_map.vertexSet());
						partner_list.sort(new Comparator<Integer>() {
							@Override
							public int compare(Integer o1, Integer o2) {
								return Integer.compare(o1.intValue(), o2.intValue());
							}
						});

						int pt = Collections.binarySearch(partner_list, indexSel);
						partner_list.remove(pt);
						Integer[] partner_arr = partner_list.toArray(new Integer[partner_list.size()]);

						int partnerBeforeTime = lookup_firstMatchFrom.get(indexSel)
								+ AbstractIndividualInterface.ONE_YEAR_INT;

						if (partnerBeforeTime > getGlobalTime()) {
							int numConnectionToAdd = (int) Math.min(dI - colSel, 
									CombinatoricsUtils.binomialCoefficient(partner_arr.length, 2));														
							
							ArrayList<String> alreadyMatched = new ArrayList<>();

							// Schedule a future partnership to schedule
							while (numConnectionToAdd > 0) {
								Integer[] selectPartner = ArrayUtilsRandomGenerator.randomSelect(partner_arr, 2,
										RNG_density);
								Arrays.sort(selectPartner);
								int sP = Collections.binarySearch(alreadyMatched, Arrays.toString(selectPartner));
								if (sP < 0) {
									alreadyMatched.add(~sP, Arrays.toString(selectPartner));
									int partnershipFormTime = getGlobalTime()
											+ RNG_density.nextInt(partnerBeforeTime - getGlobalTime());
									Integer[] newEdge = new Integer[] { selectPartner[0], selectPartner[1],
											CANDIDATE_ARRAY_SOUGHT_ANY };

									ArrayList<Integer[]> schedule_partnership_ent = schedule_partnership
											.get(partnershipFormTime);
									if (schedule_partnership_ent == null) {
										schedule_partnership_ent = new ArrayList<>();
										schedule_partnership.put(partnershipFormTime,
												schedule_partnership_ent);
									}

									int entPt = Collections.binarySearch(schedule_partnership_ent, newEdge,
											schedule_partnership_comparator);
									if (entPt < 0) {
										schedule_partnership_ent.add(~entPt, newEdge);										
									}
									numConnectionToAdd--;
								}
							}
						}
						count[dI]++;						
					}

				}

				System.out.printf("Time = %d. Unadj. density=%s. Adj. density = %s.\n", 
						getGlobalTime(), init_count, Arrays.toString(count));
			}

		}
	}

	private ContactMap contactMapInLast12Months(ContactMap baseMap) {
		if (baseMap == null) {
			return null;
		} else {
			ContactMap resMap = new ContactMap();
			Set<Integer[]> edges = baseMap.edgeSet();
			for (Integer[] e : edges) {
				if (isEdgelast12Month(e)) {
					Integer p1 = e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1];
					Integer p2 = e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2];
					if (!resMap.containsVertex(p1)) {
						resMap.addVertex(p1);
					}
					if (!resMap.containsVertex(p2)) {
						resMap.addVertex(p2);
					}
					resMap.addEdge(p1, p2, e);
				}
			}
			return resMap;
		}
	}

	private boolean hasEdgeIn12Months(ContactMap map, int pid_1, int pid_2) {
		Set<Integer[]> possibleEdges = map.getAllEdges(pid_1, pid_2);
		if (possibleEdges == null) {
			return false;
		} else {
			return possibleEdges.size() > 0;
		}

	}

	private boolean isEdgelast12Month(Integer[] e) {
		int eStart = e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME];
		int eEnd = e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]
				+ e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION];

		return eStart < getGlobalTime() && eEnd >= getGlobalTime() - AbstractIndividualInterface.ONE_YEAR_INT;
	}

}
