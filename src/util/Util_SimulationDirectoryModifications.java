/**
 * 
 */
package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.compress.archivers.sevenz.SevenZArchiveEntry;
import org.apache.commons.compress.archivers.sevenz.SevenZOutputFile;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import sim.SimulationInterface;
import sim.Simulation_ClusterModelTransmission;

/**
 * Miscellaneous helper code to modify (i.e. combine, generate) multiple
 * simulation directories
 */
public class Util_SimulationDirectoryModifications {

	public static HashMap<File, Integer> generateSeedFilesFromOptSummary(File optSummaryFile, File target_dirs_base,
			String target_dir_name_format, int numParam, int numParamSetToIncl_Max, int numParamPerSeedList_Max,
			File[] directFileCopy, String seed_list_header, boolean unqiueCombinationOnly) throws IOException {

		BufferedReader reader = new BufferedReader(new FileReader(optSummaryFile));
		String line;
		line = reader.readLine(); // Skip first line

		ArrayList<long[]> included_seed = new ArrayList<>();
		ArrayList<String[]> line_selected = new ArrayList<>();

		while ((line = reader.readLine()) != null && line_selected.size() < numParamSetToIncl_Max) {
			String[] ent = line.split(",");
			long[] seed_arr = new long[] { Long.parseLong(ent[1]), Long.parseLong(ent[2]) };
			boolean toBeInclude = true;

			if (unqiueCombinationOnly) {
				int pt = Collections.binarySearch(included_seed, seed_arr, new Comparator<long[]>() {
					@Override
					public int compare(long[] o1, long[] o2) {
						int res = 0;
						for (int i = 0; i < o1.length && res == 0; i++) {
							res = Long.compare(o1[i], o2[i]);
						}
						return res;
					}
				});
				toBeInclude = pt < 0;
				if (toBeInclude) {
					included_seed.add(~pt, seed_arr);
				}
			}
			if (toBeInclude) {
				int pt = Collections.binarySearch(line_selected, ent, new Comparator<String[]>() {
					@Override
					public int compare(String[] o1, String[] o2) {
						int res = 0;
						for (int i = 0; i < o1.length && res == 0; i++) {
							switch (i) {
							case 1: // Seed
							case 2:
								res = Long.compare(Long.parseLong(o1[i]), Long.parseLong(o2[i]));
								break;
							default:
								res = Double.compare(Double.parseDouble(o1[i]), Double.parseDouble(o2[i]));

							}
						}
						return res;
					}

				});

				if (pt < 0) {
					line_selected.add(~pt, ent);
				}
			}
		}
		reader.close();

		// Generating seed files
		final String seedFileName = "Seed_List_%d.csv";
		long last_cmap = -1;
		File lastSeedFile = null;
		int numInSeed = 0;
		int seedFile_num = 0;
		File tarDir = null;

		HashMap<File, Integer> mapping_num_seed_map = new HashMap<>();

		PrintWriter pWri_seedList = null;
		for (String[] param_set : line_selected) {
			long cMap_seed = Long.parseLong(param_set[1]);
			if (last_cmap != cMap_seed || numInSeed == numParamPerSeedList_Max) {
				if (pWri_seedList != null) {
					mapping_num_seed_map.put(lastSeedFile, numInSeed);
					pWri_seedList.close();
				}

				if (target_dir_name_format != null) {
					tarDir = tarDir == null ? new File(target_dirs_base, target_dir_name_format)
							: new File(target_dirs_base,
									String.format("%s_Extra_%d", target_dir_name_format, seedFile_num));
					tarDir.mkdirs();
				} else {
					// Group all seed list in same folder
					tarDir = target_dirs_base;
				}

				lastSeedFile = new File(tarDir, String.format(seedFileName, seedFile_num));

				pWri_seedList = new PrintWriter(lastSeedFile);
				pWri_seedList.println(seed_list_header);

				numInSeed = 0;
				seedFile_num++;
				last_cmap = cMap_seed;
			}

			for (int i = 1; i < Math.min(param_set.length, numParam + 3); i++) {
				if (i != 1) {
					pWri_seedList.print(',');
				}
				pWri_seedList.print(param_set[i]);

			}
			pWri_seedList.println();
			numInSeed++;
		}

		mapping_num_seed_map.put(lastSeedFile, numInSeed);
		pWri_seedList.close();

		// Copying other files
		for (File seedFile : mapping_num_seed_map.keySet()) {
			for (File df : directFileCopy) {
				Files.copy(df.toPath(), new File(seedFile.getParent(), df.getName()).toPath(),
						StandardCopyOption.REPLACE_EXISTING);
			}
		}

		return mapping_num_seed_map;

	}

	/**
	 * Generate multiple simulations directories based on optimisation summary file
	 * 
	 * @param optSummaryFile         Optimisation summary file.
	 * @param target_dirs_base       Target dire where the simulations directory
	 *                               will be located.
	 * @param target_dir_name_format Format of generated simulation directory.
	 * @param numParamSetToIncl      Number of parameter set to be included.
	 * @param numEntryPerSimDir      Number simulations per director
	 * @param numSimDirs             Number of directory to generate, or based on
	 *                               cMap if <0.
	 * @param prop_file_copy         Location of PROP file to copy from.
	 * @param riskGrp_def_dir        Location of directories containing risk group
	 *                               definition.
	 * @param directFileCopy         List of file to directly copy from.
	 * @param seed_list_header       Seed list header for generated Seed_List.csv
	 * 
	 * @throws FileNotFoundException
	 * @throws IOException
	 */

	public static void generateSimDirsFromOptSummary(File optSummaryFile, File target_dirs_base,
			String target_dir_name_format, int numParamSetToIncl, int numEntryPerSimDir, int numSimDirs,
			File prop_file_copy, File riskGrp_def_dir, File[] directFileCopy, String seed_list_header)
			throws FileNotFoundException, IOException {
		BufferedReader reader = new BufferedReader(new FileReader(optSummaryFile));
		String line;
		line = reader.readLine();

		HashMap<Long, ArrayList<String[]>> mapLines = new HashMap<>();

		while ((line = reader.readLine()) != null && numParamSetToIncl > 0) {
			String[] ent = line.split(",");
			Long cMap_Seed = Long.parseLong(ent[1]);
			ArrayList<String[]> linesEnt = mapLines.get(cMap_Seed);
			if (linesEnt == null) {
				linesEnt = new ArrayList<>();
				mapLines.put(cMap_Seed, linesEnt);
			}
			linesEnt.add(Arrays.copyOfRange(ent, 1, ent.length - 1));
			numParamSetToIncl--;
		}
		reader.close();

		ArrayList<ArrayList<String[]>> targetDir_collection = new ArrayList<>();

		if (numSimDirs > 0) {

			numEntryPerSimDir = numParamSetToIncl / numSimDirs;
			while (numEntryPerSimDir * numSimDirs < numParamSetToIncl) {
				numEntryPerSimDir++;
			}

			ArrayList<ArrayList<String[]>> entry_collection = new ArrayList<>(mapLines.values());
			Comparator<ArrayList<String[]>> greedy_Comparator = new Comparator<ArrayList<String[]>>() {
				@Override
				public int compare(ArrayList<String[]> o1, ArrayList<String[]> o2) {
					return -Integer.compare(o1.size(), o2.size());
				}
			};
			Collections.sort(entry_collection, greedy_Comparator);

			for (int i = 0; i < numSimDirs; i++) {
				targetDir_collection.add(new ArrayList<>());
			}

			ArrayList<String[]> largest_entry = entry_collection.get(0);
			while (largest_entry.size() > 0) {
				ArrayList<String[]> targetDir_ent = targetDir_collection.get(targetDir_collection.size() - 1);
				for (int i = targetDir_ent.size(); i < numEntryPerSimDir && largest_entry.size() > 0; i++) {
					targetDir_ent.add(largest_entry.remove(0));
				}
				Collections.sort(entry_collection, greedy_Comparator);
				Collections.sort(targetDir_collection, greedy_Comparator);
				largest_entry = entry_collection.get(0);
			}
		} else {
			for (Long cMapSeed : mapLines.keySet()) {
				ArrayList<String[]> tar_ent = new ArrayList<String[]>(numEntryPerSimDir);
				targetDir_collection.add(tar_ent);
				ArrayList<String[]> param_ent = mapLines.get(cMapSeed);
				for (String[] param_str : param_ent) {
					if (tar_ent.size() >= numEntryPerSimDir) {
						tar_ent = new ArrayList<String[]>(numEntryPerSimDir);
						targetDir_collection.add(tar_ent);
					}
					tar_ent.add(param_str);
				}
			}
			numSimDirs = targetDir_collection.size();
		}

		File[] tarDir = new File[numSimDirs];

		for (int f = 0; f < tarDir.length; f++) {
			tarDir[f] = new File(target_dirs_base, String.format(target_dir_name_format, f));
			tarDir[f].mkdirs();
			PrintWriter pWri = new PrintWriter(new File(tarDir[f], "Seed_List.csv"));
			pWri.println(seed_list_header);

			ArrayList<String[]> linesEnt = targetDir_collection.get(f);
			ArrayList<Long> cMapSeedArr = new ArrayList<>();
			if (linesEnt != null) {
				for (String[] ent : linesEnt) {
					long cMapSeed = Long.parseLong(ent[0]);
					int k = Collections.binarySearch(cMapSeedArr, cMapSeed);
					if (k < 0) {
						cMapSeedArr.add(~k, cMapSeed);
					}
					boolean newLine = true;
					for (String s : ent) {
						if (!newLine) {
							pWri.print(",");
						}
						pWri.print(s);
						newLine = false;
					}
					pWri.println();
				}
			}
			pWri.close();

			for (long cMapSeed : cMapSeedArr) {
				File riskGrp_File = new File(riskGrp_def_dir,
						String.format(Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP, cMapSeed));

				if (riskGrp_File.exists()) {
					Files.copy(riskGrp_File.toPath(), new File(tarDir[f], riskGrp_File.getName()).toPath(),
							StandardCopyOption.REPLACE_EXISTING);
				}
			}

			// Copying prop file

			BufferedReader prop_reader = new BufferedReader(new FileReader(prop_file_copy));
			PrintWriter prop_writer = new PrintWriter(new File(tarDir[f], prop_file_copy.getName()));

			String propLine;
			while ((propLine = prop_reader.readLine()) != null) {
				if (propLine.startsWith("<entry key=\"PROP_NUM_SIM_PER_SET\">")) {
					propLine = String.format("<entry key=\"PROP_NUM_SIM_PER_SET\">%d</entry>", linesEnt.size());
				} else if (propLine.startsWith("<entry key=\"PROP_USE_PARALLEL\">")) {
					propLine = String.format("<entry key=\"PROP_USE_PARALLEL\">%d</entry>", linesEnt.size());
				}
				prop_writer.println(propLine);
			}

			prop_reader.close();
			prop_writer.close();

			for (File df : directFileCopy) {
				Files.copy(df.toPath(), new File(tarDir[f], df.getName()).toPath(),
						StandardCopyOption.REPLACE_EXISTING);
			}

			System.out.printf("In %s -> %d entries with cMap seed of %s\n", tarDir[f].getName(), linesEnt.size(),
					Arrays.toString(cMapSeedArr.toArray(new Long[cMapSeedArr.size()])));

		}

		System.out.printf("Splited seed file generated at %s.\n", target_dirs_base.getAbsolutePath());
	}

	public static void generateUniqueZip(File tarDirPath, String zip_file_name, Pattern match_pattern)
			throws IOException {

		File[] targetDirs = tarDirPath.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory();
			}
		});

		ArrayList<String> includedFilename = new ArrayList<>();
		ArrayList<File> toBeRemoved = new ArrayList<>();
		ArrayList<File> toBeAdded = new ArrayList<>();

		for (File tarDir : targetDirs) {
			File[] zip_candidates = tarDir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return match_pattern.matcher(pathname.getName()).matches();
				}
			});
			for (File zipCandidate : zip_candidates) {
				int pt = Collections.binarySearch(includedFilename, zipCandidate.getName());
				if (pt < 0) {
					includedFilename.add(~pt, zipCandidate.getName());

					toBeAdded.add(zipCandidate);

				}
				toBeRemoved.add(zipCandidate);
			}
		}
		File zip_tar_file = new File(tarDirPath, zip_file_name);

		if (toBeAdded.size() > 0) {
			SevenZOutputFile outputZip = new SevenZOutputFile(zip_tar_file);
			SevenZArchiveEntry entry;
			FileInputStream fIn;
			for (File addFile : toBeAdded) {
				entry = outputZip.createArchiveEntry(addFile, addFile.getName());
				outputZip.putArchiveEntry(entry);
				fIn = new FileInputStream(addFile);
				outputZip.write(fIn);
				outputZip.closeArchiveEntry();
				fIn.close();
			}
			outputZip.close();
		}

		for (File rmFile : toBeRemoved) {
			Files.delete(rmFile.toPath());
		}
		System.out.printf("Generating unique zip for %s completed. %d files removed, with %d files add to %s.\n",
				tarDirPath.getAbsolutePath(),
				toBeRemoved.size(),
				includedFilename.size(), zip_tar_file.getName());

	}

	/**
	 * Combine simulation outputs directory
	 * 
	 * @param tarDirPath Directory where simulation output is stored.
	 * @param extraDir   Directory where extra directories are stored
	 * @throws IOException
	 * @throws FileNotFoundException
	 * 
	 */
	public static void combineSimOutput(File tarDirPath, File extraDir) throws IOException, FileNotFoundException {

		String[] ignoreAttr = new String[] { "PROP_CONTACT_MAP_LOC", "PROP_NUM_SIM_PER_SET", "PROP_USE_PARALLEL" };
		String[][] post_zip_patterns = new String[][] { new String[] {
				Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP.replaceAll("_%d", "") + ".7z",
				Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP.replaceAll("%d",
						"(-{0,1}\\\\d+)"), },
				new String[] { "Seed_List.csv.zip", "Seed_List(_\\d+).csv", },
				new String[] { "simSpecificSwitch.prop.zip", "simSpecificSwitch(_\\d+).prop", },
				new String[] { "dx_Bali.prop.zip", "dx_Bali(_\\d+).prop", }, };

		Arrays.sort(ignoreAttr);

		File[] targetDirs = tarDirPath.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory() && !pathname.equals(extraDir);
			}
		});

		Pattern resPattern = Pattern.compile("\\[(.*)\\](.*)_(-{0,1}\\d+)_(-{0,1}\\d+).csv");

		for (File tarDir : targetDirs) {

			File[] srcDirs = extraDir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pathname.isDirectory() && pathname.getName().startsWith(tarDir.getName());
				}
			});

			System.out.printf("Combining %d directories to %s\n", srcDirs.length, tarDir.getAbsolutePath());

			for (int f = 0; f < srcDirs.length; f++) {
				File srcDir = srcDirs[f];
				File[] srcEnts = srcDir.listFiles();
				for (File srcEnt : srcEnts) {
					File tarEnt = new File(tarDir, srcEnt.getName());

					if (tarEnt.exists()) {
						if (tarEnt.getName().equals(SimulationInterface.FILENAME_PROP)) {
							try {
								Document xml_doc = PropValUtils.parseXMLFile(tarEnt);
								NodeList nList_entry = xml_doc.getElementsByTagName("entry");
								HashMap<String, String> entryMap = new HashMap<>();

								for (int entId = 0; entId < nList_entry.getLength(); entId++) {
									Element entryElement = (Element) nList_entry.item(entId);
									entryMap.put(entryElement.getAttribute("key"),
											entryElement.getTextContent().trim());
								}

								xml_doc = PropValUtils.parseXMLFile(srcEnt);
								nList_entry = xml_doc.getElementsByTagName("entry");

								for (int entId = 0; entId < nList_entry.getLength(); entId++) {
									Element entryElement = (Element) nList_entry.item(entId);
									String key = entryElement.getAttribute("key");
									if (Arrays.binarySearch(ignoreAttr, key) < 0) {
										if (!entryElement.getTextContent().trim().equals(entryMap.get(key))) {
											System.out.printf("Warning: Entry with attr = \"%s\" does not match.\n",
													key);
											System.out.printf(" Tar = \"%s\"\n", entryMap.get(key));
											System.out.printf(" Src = \"%s\"\n", entryElement.getTextContent().trim());
										}
									}
								}

							} catch (ParserConfigurationException | SAXException ex) {
								System.out.printf("%s and/or %s are invalid XML file - file not compared.\n",
										srcEnt.getAbsolutePath(), tarEnt.getAbsolutePath());

							}
						} else {
							String srcName = srcEnt.getName();
							if (srcName.equals(tarEnt.getName())) {
								int first_dot = srcName.indexOf('.');
								if (first_dot < 0) {
									first_dot = srcName.length();
								}
								tarEnt = new File(tarDir, String.format("%s_%d%s", srcName.substring(0, first_dot), f,
										srcName.substring(first_dot, srcName.length())));
							}
						}
					}

					if (!srcEnt.getName().equals(SimulationInterface.FILENAME_PROP)) {
						Files.move(srcEnt.toPath(), tarEnt.toPath());
					}
				}

			}

			// Clean up duplicate risk map
			Pattern riskGrp_pattern = Pattern.compile(Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP
					.replaceAll("%d", "(-{0,1}\\\\d+)"));
			File[] riskGrp_files = tarDir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return riskGrp_pattern.matcher(pathname.getName()).matches();
				}
			});
			for (File riskGrp_file : riskGrp_files) {
				String prefix = riskGrp_file.getName().substring(0, riskGrp_file.getName().indexOf('.')) + "_";

				File[] dup_riskGrp_files = tarDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return pathname.getName().startsWith(prefix);
					}
				});

				if (dup_riskGrp_files.length > 0) {
					System.out.printf("Cleaning up Risk Group maping for %s.\n", riskGrp_file.getName());

					for (File dup_riskGrp_file : dup_riskGrp_files) {
						Files.delete(dup_riskGrp_file.toPath());
					}
				}
			}

			// Zip up CSV
			for (String[] zip_patten : post_zip_patterns) {
				String zip_file_name = zip_patten[0];
				Pattern match_pattern = Pattern.compile(zip_patten[1]);

				File[] zip_candidates = tarDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return match_pattern.matcher(pathname.getName()).matches();
					}
				});

				if (zip_candidates.length > 0) {
					File zip_tar_file = new File(tarDir, zip_file_name);
					SevenZOutputFile outputZip = new SevenZOutputFile(zip_tar_file);

					SevenZArchiveEntry entry;
					FileInputStream fIn;

					for (int fI = 0; fI < zip_candidates.length; fI++) {
						entry = outputZip.createArchiveEntry(zip_candidates[fI], zip_candidates[fI].getName());
						outputZip.putArchiveEntry(entry);
						fIn = new FileInputStream(zip_candidates[fI]);
						outputZip.write(fIn);
						outputZip.closeArchiveEntry();
						fIn.close();
					}
					outputZip.close();

					for (int fI = 0; fI < zip_candidates.length; fI++) {
						Files.delete(zip_candidates[fI].toPath());
					}
					System.out.printf("%d files zipped to %s.\n", zip_candidates.length, zip_tar_file.getName());
				}
			}

			// Check if original CSV exist and zip it up
			File[] resCSV = tarDir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return resPattern.matcher(pathname.getName()).matches();
				}
			});

			while (resCSV.length != 0) {
				Matcher m = resPattern.matcher(resCSV[0].getName());
				m.matches();

				Pattern subCSVPatten = Pattern
						.compile(String.format("\\[%s\\]%s_%s_(-{0,1}\\d+).csv", m.group(1), m.group(2), m.group(3)));
				File subZip = new File(tarDir, String.format("%s_%s.csv.7z", m.group(2), m.group(3)));

				if (subZip.exists()) {
					System.out.printf("Error: Target zip %s already exist. Target zip are to be renamed.\n",
							subZip.getName());
					int extra = 0;
					while (subZip.exists()) {
						subZip = new File(tarDir, String.format("%s_%s_%d.csv.7z", m.group(2), m.group(3), extra));
						extra++;
					}

				}

				if (!subZip.exists()) {
					System.out.printf("Generating %s ....", subZip.getAbsolutePath());
					Simulation_ClusterModelTransmission.zipSelectedOutputs(tarDir, subZip.getName(), subCSVPatten,
							true);
					System.out.println("Done");
				}

				resCSV = tarDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return resPattern.matcher(pathname.getName()).matches();
					}
				});
			}

		}
	}
	
	/**
	 * Extract lines from seed files, either from CSV or zip 
	 * @param seedFiles
	 * @return 
	 * @throws IOException
	 * @throws FileNotFoundException
	 */

	public static ArrayList<String[]> extractSeedFilesEntries(File[] seedFiles)
			throws IOException, FileNotFoundException {
		ArrayList<String[]> seedFileTxtArray = new ArrayList<>();
	
		for (File seedFile : seedFiles) {
			if (seedFile.exists()) {
				if (seedFile.getName().endsWith(".zip")) {
					HashMap<String, ArrayList<String[]>> entMap = util.Util_7Z_CSV_Entry_Extract_Callable
							.extractedLinesFrom7Zip(seedFile);
					for (Entry<String, ArrayList<String[]>> ent : entMap.entrySet()) {
						String[] lines;
						ArrayList<String[]> src = ent.getValue();
						lines = new String[src.size()];
						for (int lineNum = 0; lineNum < lines.length; lineNum++) {
							StringBuilder line_builder = new StringBuilder();
							String[] lineEnt = src.get(lineNum);
							for (int j = 0; j < lineEnt.length; j++) {
								if (j != 0) {
									line_builder.append(',');
								}
								line_builder.append(lineEnt[j]);
							}
							lines[lineNum] = line_builder.toString();
						}
						seedFileTxtArray.add(lines);
					}
				} else {
					String[] lines;
					lines = util.Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(seedFile);
					seedFileTxtArray.add(lines);
				}
			}
		}
		return seedFileTxtArray;
	}

}
