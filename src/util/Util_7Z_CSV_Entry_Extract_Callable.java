package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.StringReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.compress.archivers.sevenz.SevenZArchiveEntry;
import org.apache.commons.compress.archivers.sevenz.SevenZFile;

public class Util_7Z_CSV_Entry_Extract_Callable implements Callable<Map<String, long[]>> {

	final private File BaseDir;
	final private String[] CSV_Filenames_Pattern;
	final private String CMapId;
	final private int[][] Entry_Pairs; // {[File_id, colNum}

	public Util_7Z_CSV_Entry_Extract_Callable(File baseDir, String cMapId, String[] csv_Filenames_Pattern,
			int[][] entry_Pairs) {
		super();
		BaseDir = baseDir;
		CMapId = cMapId;
		CSV_Filenames_Pattern = csv_Filenames_Pattern;
		Entry_Pairs = entry_Pairs;
	}

	@Override
	public Map<String, long[]> call() throws Exception {
		final int BUFFER = 2048;
		final byte[] buf = new byte[BUFFER];
		HashMap<String, long[]> res = new HashMap<>(); // key = cMapSeed_SimSeed, value =

		for (int file_id = 0; file_id < CSV_Filenames_Pattern.length; file_id++) {
			Pattern pattern_csv_file = Pattern
					.compile(CSV_Filenames_Pattern[file_id].replaceAll("%d", "(-{0,1}(?!0)\\\\d+)"));

			File archiveFile = new File(BaseDir,
					(CSV_Filenames_Pattern[file_id].replaceFirst("_%d", "") + ".7z").replaceAll("%d", CMapId));

			SevenZFile archive7Z = new SevenZFile(archiveFile);
			SevenZArchiveEntry ent;
			while ((ent = archive7Z.getNextEntry()) != null) {
				Matcher m = pattern_csv_file.matcher(ent.getName());
				if (m.matches()) {
					String res_key = String.format("%s_%s", m.group(1), m.group(2));
					long[] res_ent = res.get(res_key);
					if (res_ent == null) {
						res_ent = new long[Entry_Pairs.length];
						Arrays.fill(res_ent, -1);
						res.put(res_key, res_ent);
					}
					StringBuilder txt_entries = new StringBuilder();
					int count;
					while ((count = archive7Z.read(buf, 0, BUFFER)) != -1) {
						txt_entries.append(new String(Arrays.copyOf(buf, count)));
					}
					String line, lastline = "";
					BufferedReader lines = new BufferedReader(new StringReader(txt_entries.toString()));
					while ((line = lines.readLine()) != null) {
						if (line.length() > 0) {
							lastline = line;
						}
					}
					lines.close();
					String[] lastlineEnt = lastline.split(",");
					for (int k = 0; k < Entry_Pairs.length; k++) {
						if (Entry_Pairs[k][0] == file_id) {							
							try {
								res_ent[k] = Long.parseLong(lastlineEnt[Entry_Pairs[k][1]]);
							} catch(ArrayIndexOutOfBoundsException ex) {
								System.err.printf("Error in last line in %s : %s\n",
										BaseDir.getName(), ent.getName());
								if(res_ent[k] != -1) {									
									res_ent[k] = -1;
								}
								
							}
						}
					}
				}
			}
			archive7Z.close();

		}
		return res;
	}

}