package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import person.AbstractIndividualInterface;
import relationship.ContactMap;
import sim.Abstract_Runnable_ClusterModel;
import sim.Runnable_ClusterModel_ContactMap_Generation_MultiMap;
import sim.Simulation_ClusterModelGeneration;
import sim.Simulation_ClusterModelTransmission;

public class Util_ContactMap_Adjust {

	File dirSrcBase, dirSrcNew, dirRes;

	final Pattern PATTERN_CONTACT_MAP = Pattern.compile(
			Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)"));
	
	
	public Util_ContactMap_Adjust(String dirSrcBase, String dirSrcNew, String dirRes) {
		super();
		this.dirSrcBase = new File(dirSrcBase);
		this.dirSrcNew = new File(dirSrcNew);
		this.dirRes = new File(dirRes);
	}

	

	public Util_ContactMap_Adjust(File dirSrcBase, File dirSrcNew, File dirRes) {
		super();
		this.dirSrcBase = dirSrcBase;
		this.dirSrcNew = dirSrcNew;
		this.dirRes = dirRes;
	}

	public void combineMaps() throws FileNotFoundException, IOException {
		File[] base_map_arr = dirSrcBase.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return PATTERN_CONTACT_MAP.matcher(pathname.getName()).matches();
			}
		});

		for (File base_map : base_map_arr) {

			Matcher m = PATTERN_CONTACT_MAP.matcher(base_map.getName());
			m.matches();
			long mapSeed = Long.parseLong(m.group(1));

			File[] new_map_arr = dirSrcNew.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return Pattern.matches(Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP
							.replaceFirst("%d", Long.toString(mapSeed)).replaceFirst("%d", "(-{0,1}(?!0)\\\\d+)"),
							pathname.getName());
				}
			});

			if (new_map_arr.length > 0) {

				System.out.printf("Matching %d file(s) with contact map located in %s.\n", new_map_arr.length,
						base_map.getName());

				BufferedReader reader;
				StringWriter lines;
				PrintWriter pWri;
				String line;

				int maxEdgeStart = -1;

				reader = new BufferedReader(new FileReader(base_map));
				lines = new StringWriter();
				pWri = new PrintWriter(lines);

				while ((line = reader.readLine()) != null) {
					String[] ent = line.split(",");
					maxEdgeStart = Math.max(maxEdgeStart,
							Integer.parseInt(ent[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]));
					pWri.println(line);
				}
				reader.close();
				pWri.close();

				System.out.printf("Last edge in base map starts at %d.\n", maxEdgeStart);

				ContactMap base = ContactMap.ContactMapFromFullString(lines.toString());

				lines.close();

				for (File new_map : new_map_arr) {
					reader = new BufferedReader(new FileReader(new_map));
					while ((line = reader.readLine()) != null) {
						String[] ent = line.split(",");
						Integer[] newEdge = new Integer[ent.length];

						for (int i = 0; i < newEdge.length; i++) {
							newEdge[i] = Integer.parseInt(ent[i]);
						}

						Integer[] targetEdge = null;
						Set<Integer[]> edgeSet = base.getAllEdges(
								newEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1],
								newEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2]);

						if (edgeSet != null) {

							for (Integer[] testEdge : edgeSet) {
								if (testEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]
										.equals(newEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME])) {
									targetEdge = testEdge;
								}
							}
						}

						if (targetEdge != null) {
							targetEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] = Math.max(
									targetEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION],
									newEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION]);
						} else {
							maxEdgeStart = Math.max(maxEdgeStart,
									newEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]);
							targetEdge = newEdge;
							
							if(!base.containsVertex(targetEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1])) {
								base.addVertex(newEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1]);
							}
							
							if(!base.containsVertex(targetEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2])) {
								base.addVertex(newEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2]);
							}
							
							
							base.addEdge(targetEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1],
									targetEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2], targetEdge);

						}

					}
				}

				File targetFile = new File(dirRes, Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP
						.replaceFirst("%d", Long.toString(mapSeed)).replaceFirst("%d", Integer.toString(maxEdgeStart)));

				pWri = new PrintWriter(targetFile);
				pWri.println(base.toFullString());
				pWri.close();

				System.out.printf("Combine contact map %s generated.\n", targetFile.getName());

			}

		}

	}



	public static void convertContactMapToMultiMap(File contactMapDir, File tarDir)
			throws FileNotFoundException, IOException {
		final Pattern pattern_all_ContactMap = Pattern.compile(
				Simulation_ClusterModelTransmission.FILENAME_FORMAT_ALL_CMAP.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)"));
			
		tarDir.mkdirs();	
	
		File[] preGenClusterMap = contactMapDir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isFile() && pattern_all_ContactMap.matcher(pathname.getName()).matches();
	
			}
		});
	
		System.out.printf("Converting %d contact map found in %s...\n", preGenClusterMap.length,
				contactMapDir.getAbsolutePath());
	
		PrintWriter pwri;
		for (File file_cmap : preGenClusterMap) {
			Matcher m = pattern_all_ContactMap.matcher(file_cmap.getName());
			if (m.find()) {
				ArrayList<Integer> pid_list = new ArrayList<Integer>();
				int minTime = Integer.MAX_VALUE;
				int maxTime = 0;
				int default_enter_age = 18 * AbstractIndividualInterface.ONE_YEAR_INT;
				
				long cMap_seed = Long.parseLong(m.group(1));
				HashMap<Integer, ArrayList<int[]>> futureMap = new HashMap<>();
				BufferedReader cmap_reader = new BufferedReader(new FileReader(file_cmap));												
				
				String cmap_line;
			
				File convert_cmap_file = new File(tarDir, String.format(Runnable_ClusterModel_ContactMap_Generation_MultiMap.MAPFILE_FORMAT, 0, cMap_seed));
				pwri = new PrintWriter(convert_cmap_file);
	
				while ((cmap_line = cmap_reader.readLine()) != null) {
					String[] cmap_entry = cmap_line.split(",");
	
					int pid_1 = Integer.parseInt(cmap_entry[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1]);
					int pid_2 = Integer.parseInt(cmap_entry[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2]);
					int pre_e_start = -1;
	
					for (int pid : new int[] { pid_1, pid_2 }) {
						int pt = Collections.binarySearch(pid_list, pid);
						if (pt < 0) {
							pid_list.add(~pt, pid);
						}
					}
	
					for (int i = Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME; i < cmap_entry.length; i += 2) {
						int e_start = Integer.parseInt(cmap_entry[i]);
						int e_length = Integer.parseInt(cmap_entry[i + 1]);
						
						minTime = Math.min(e_start, minTime);
						maxTime = Math.max(e_start+e_length, maxTime);
	
						if (pre_e_start != e_start) {
							Integer[] storedEdgeTimes = futureMap.keySet().toArray(new Integer[0]);
							if (storedEdgeTimes.length > 0) {
								Arrays.sort(storedEdgeTimes);
								int pt = 0;
								while (pt < storedEdgeTimes.length 
										&& storedEdgeTimes[pt] <= e_start) {
									ArrayList<int[]> storedEdges = futureMap.remove(storedEdgeTimes[pt]);
									for (int[] edges : storedEdges) {
										pwri.printf("%d,%d,%d,%d\n", edges[0], edges[1],edges[2],edges[3]);
									}
									pt++;
								}
	
							}
							pre_e_start = e_start;
						}
						if (i != Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME) {
							ArrayList<int[]> ent = futureMap.get(e_start);
							if (ent == null) {
								ent = new ArrayList<>();
								futureMap.put(e_start, ent);
							}
							ent.add(new int[] { pid_1, pid_2, e_start, e_length });
						}else {
							pwri.printf("%d,%d,%d,%d\n", pid_1, pid_2, e_start, e_length);
						}
					}
				}
				
				pwri.close();
				cmap_reader.close();
				
				System.out.printf("Convert contact map to %s\n", convert_cmap_file.getAbsolutePath());
	
				// Print pop stat
				
				File pop_stat_file = new File(tarDir,
						String.format(Runnable_ClusterModel_ContactMap_Generation_MultiMap.POPSTAT_FORMAT, cMap_seed));
				pwri = new PrintWriter(pop_stat_file);
				pwri.println("ID,GRP,ENTER_POP_AGE,ENTER_POP_AT,EXIT_POP_AT,HAS_REG_PARTNER_UNTIL");
	
				for (Integer pid : pid_list) {
					pwri.printf("%d,0,%d,%d,%d,-1\n", pid, default_enter_age, minTime, maxTime);
				}
	
				pwri.close();
				
				System.out.printf("Generated default pop stat file at %s \n", pop_stat_file.getAbsolutePath());
	
			}
	
		}
	}

}
