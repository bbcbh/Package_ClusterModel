package sim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Set;
import java.util.concurrent.Callable;

import relationship.ContactMap;
import relationship.TransmissionMap;

public abstract class Abstract_Runnable_ClusterModel implements Runnable {

	protected String runnableId = null;
	protected File baseDir = null;
	
	// Population Index
	// Offset by ID under POP_STAT_%d.csv
	public static final int POP_INDEX_GRP = 0;	
	public static final int POP_INDEX_ENTER_POP_AGE = POP_INDEX_GRP + 1;
	public static final int POP_INDEX_ENTER_POP_AT = POP_INDEX_ENTER_POP_AGE + 1;
	public static final int POP_INDEX_EXIT_POP_AT = POP_INDEX_ENTER_POP_AT + 1;
	public static final int POP_INDEX_HAS_REG_PARTNER_UNTIL = POP_INDEX_EXIT_POP_AT + 1;
	public static final int LENGTH_POP_ENTRIES = POP_INDEX_HAS_REG_PARTNER_UNTIL + 1;

	// Edge format in contact map: {p1, p2, start time_1, duration_1, ...}
	// CONTACT_MAP_EDGE_DURATION = -1 -> Removal
	// CONTACT_MAP_EDGE_DURATION = 0 -> Update
	public static final int CONTACT_MAP_EDGE_P1 = 0;
	public static final int CONTACT_MAP_EDGE_P2 = CONTACT_MAP_EDGE_P1 + 1;
	public static final int CONTACT_MAP_EDGE_START_TIME = CONTACT_MAP_EDGE_P2 + 1;
	public static final int CONTACT_MAP_EDGE_DURATION = CONTACT_MAP_EDGE_START_TIME + 1;
	public static final int LENGTH_CONTACT_MAP_EDGE = CONTACT_MAP_EDGE_DURATION + 1;

	// Edge format in transmission map
	public static final int TRANS_MAP_EDGE_INFECTIOUS = 0;
	public static final int TRANS_MAP_EDGE_SUSCEPTIBLE = TRANS_MAP_EDGE_INFECTIOUS + 1;
	public static final int TRANS_MAP_EDGE_START_TIME = TRANS_MAP_EDGE_SUSCEPTIBLE + 1;
	public static final int TRANS_MAP_EDGE_ACT_INVOLVED = TRANS_MAP_EDGE_START_TIME + 1;
	public static final int LENGTH_TRANS_MAP_EDGE = TRANS_MAP_EDGE_ACT_INVOLVED + 1;

	public String getRunnableId() {
		return runnableId;
	}

	public void setRunnableId(String runnableId) {
		this.runnableId = runnableId;
	}

	public void setBaseDir(File baseDir) {
		this.baseDir = baseDir;
	}

	public abstract Object[] getRunnable_fields();	


	public static ContactMap generateContactMapAcrossTimeRange(ContactMap baseMap, int[] time_range) {
		if (baseMap == null) {
			return null;
		} else {
			ContactMap resMap = new ContactMap();
			Set<Integer[]> edges = baseMap.edgeSet();
			for (Integer[] e : edges) {
				if (isEdgeBetweenTimeRange(e, time_range)) {
					Integer p1 = e[CONTACT_MAP_EDGE_P1];
					Integer p2 = e[CONTACT_MAP_EDGE_P2];
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

	public static boolean isEdgeBetweenTimeRange(Integer[] e, int[] time_range) {
		int eStart = e[CONTACT_MAP_EDGE_START_TIME];
		int eEnd = e[CONTACT_MAP_EDGE_START_TIME] + e[CONTACT_MAP_EDGE_DURATION];
		return eStart < time_range[1] && eEnd >= time_range[0];
	}

	public static final Callable<ContactMap> generateContactMapCallable(File cMap_file) {
		Callable<ContactMap> callable = new Callable<ContactMap>() {
			@Override
			public ContactMap call() throws Exception {
				long tic = System.currentTimeMillis();
				StringWriter cMap_str = new StringWriter();
				PrintWriter pWri = new PrintWriter(cMap_str);
				BufferedReader reader = new BufferedReader(new FileReader(cMap_file));
				String line;
				while ((line = reader.readLine()) != null) {
					pWri.println(line);
				}

				pWri.close();
				reader.close();
				ContactMap map = ContactMap.ContactMapFromFullString(cMap_str.toString());
				System.out.printf("Contact map read from %s. Time required = %.3fs\n", cMap_file.getName(),
						(System.currentTimeMillis() - tic) / 1000f);

				return map;
			}

		};
		return callable;
	}

	public static final Callable<ContactMap> generateContactMapCallable(String cMap_str) {
		Callable<ContactMap> callable = new Callable<ContactMap>() {
			@Override
			public ContactMap call() throws Exception {
				return ContactMap.ContactMapFromFullString(cMap_str);
			}
		};
		return callable;
	}

	public static final Callable<ArrayList<Integer[]>> generateMapEdgeArray(TransmissionMap tMap) {
		Callable<ArrayList<Integer[]>> callable = new Callable<ArrayList<Integer[]>>() {
			@Override
			public ArrayList<Integer[]> call() throws Exception {
				Comparator<Integer[]> edge_cmp = generateMapEdgeComparator(INDEX_ORDER_TRANSMISSION_MAP);
				ArrayList<Integer[]> edges = new ArrayList<>(tMap.edgeSet().size());
				for (Integer[] ent : tMap.edgeSet()) {
					int key = Collections.binarySearch(edges, ent, edge_cmp);
					if (key >= 0) {
						System.err.printf("Warning from Transmission Map #%d: Edge %s already in list. Edge skipped.\n",
								tMap.getId(), Arrays.toString(ent));
					} else {
						edges.add(~key, ent);
					}

				}

				return edges;
			}

		};
		return callable;
	}

	public static final Callable<ArrayList<Integer[]>> generateMapEdgeArray(ContactMap cMap) {

		Callable<ArrayList<Integer[]>> callable = new Callable<ArrayList<Integer[]>>() {
			@Override
			public ArrayList<Integer[]> call() throws Exception {
				Comparator<Integer[]> edge_cmp = generateMapEdgeComparator(INDEX_ORDER_CONTACT_MAP);
				ArrayList<Integer[]> edges = new ArrayList<>(cMap.edgeSet().size());
				for (Integer[] ent : cMap.edgeSet()) {
					int startTime = Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME;
					int duration = Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION;
					while (duration < ent.length) {
						Integer[] e = new Integer[] { ent[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1],
								ent[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2], ent[startTime],
								ent[duration] };
						int key = Collections.binarySearch(edges, e, edge_cmp);

						if (key >= 0) {
							// System.err.printf("Warning from ContactMap #%d: Edge %s from %s already in
							// list. Edge skipped.\n",
							// cMap.getId(), Arrays.toString(e), Arrays.toString(ent));
						} else {
							edges.add(~key, e);
						}
						startTime += 2;
						duration += 2;
					}
				}
				return edges;
			}
		};
		return callable;
	}

	@Deprecated
	public static final Callable<ArrayList<Integer[]>> generateContactEdgeArray(File cMap_file) {
		Callable<ArrayList<Integer[]>> callable = new Callable<ArrayList<Integer[]>>() {
			@Override
			public ArrayList<Integer[]> call() throws Exception {
				StringWriter cMap_str = new StringWriter();
				PrintWriter pWri = new PrintWriter(cMap_str);
				BufferedReader reader = new BufferedReader(new FileReader(cMap_file));
				String line;
				int minNumLine = 0;
				while ((line = reader.readLine()) != null) {
					pWri.println(line);
					minNumLine++;
				}
				pWri.close();
				reader.close();
				return generateContactEdgeArray(cMap_str.toString(), minNumLine).call();
			}
		};
		return callable;
	}

	@Deprecated
	public static final Callable<ArrayList<Integer[]>> generateContactEdgeArray(String cMap_str, int minSize) {
		Callable<ArrayList<Integer[]>> callable = new Callable<ArrayList<Integer[]>>() {

			@Override
			public ArrayList<Integer[]> call() throws Exception {
				ArrayList<Integer[]> edges = new ArrayList<>(minSize);
				BufferedReader lines = new BufferedReader(new StringReader(cMap_str));
				String line;
				Comparator<Integer[]> edge_cmp = generateMapEdgeComparator(INDEX_ORDER_CONTACT_MAP);

				while ((line = lines.readLine()) != null) {
					String[] ent = line.split(",");
					Integer p1 = Integer.parseInt(ent[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1]);
					Integer p2 = Integer.parseInt(ent[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2]);
					int startTime = Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME;
					int duration = Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION;

					while (duration < ent.length) {
						Integer[] e = new Integer[] { p1, p2, Integer.parseInt(ent[startTime]),
								Integer.parseInt(ent[duration]), };

						int key = Collections.binarySearch(edges, e, edge_cmp);
						edges.add(~key, e);

						startTime += 2;
						duration += 2;
					}

				}

				return edges;
			}

		};

		return callable;
	}

	private final static int[] INDEX_ORDER_CONTACT_MAP = new int[] {
			Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME,
			Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION,
			Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1, Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2, };

	private final static int[] INDEX_ORDER_TRANSMISSION_MAP = new int[] {
			Abstract_Runnable_ClusterModel.TRANS_MAP_EDGE_START_TIME,
			Abstract_Runnable_ClusterModel.TRANS_MAP_EDGE_ACT_INVOLVED,
			Abstract_Runnable_ClusterModel.TRANS_MAP_EDGE_INFECTIOUS,
			Abstract_Runnable_ClusterModel.TRANS_MAP_EDGE_SUSCEPTIBLE };

	private static Comparator<Integer[]> generateMapEdgeComparator(final int[] CMP_INDEX_ORDER) {
		Comparator<Integer[]> edge_cmp = new Comparator<Integer[]>() {
			@Override
			public int compare(Integer[] o1, Integer[] o2) {
				int cmp = 0;
				for (int cmpIndex : CMP_INDEX_ORDER) {
					cmp = Integer.compare(o1[cmpIndex], o2[cmpIndex]);
					if (cmp != 0) {
						break;
					}
				}
				return cmp;
			}
		};
		return edge_cmp;
	}

}
