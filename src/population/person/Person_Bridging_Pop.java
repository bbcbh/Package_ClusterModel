package population.person;

import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import infection.AbstractInfection;
import person.AbstractIndividualInterface;

public class Person_Bridging_Pop implements AbstractIndividualInterface {

	private Object[] fields = new Object[] {
			// FIELD_ID
			-1,
			// FIELD_GENDER
			0,
			// FIELD_AGE
			18.0 * AbstractIndividualInterface.ONE_YEAR_INT,
			// FIELD_ENTER_POP_AT_AGE
			18.0 * AbstractIndividualInterface.ONE_YEAR_INT,
			// FIELD_ENTER_POP_AT_TIME
			0,
			// FIELD_MAX_REGULAR_PARTNER
			0,
			// FIELD_MAX_CASUAL_PARTNERS_6_MONTHS
			0,
			// FIELD_INFECT_STAT
			new int[0],
			// FIELD_LAST_INFECTED_AT_AGE
			new int[0],
			// FIELD_TIME_UNTIL_NEXT_STAGE
			new int[0],
			// FIELD_NUM_PARTNER
			new int[2],
			// FIELD_LAST_ACT_INFECTIOUS
			new boolean[0],

	};

	// Int field
	public static final int FIELD_ID = 0;
	public static final int FIELD_GENDER = FIELD_ID + 1;
	public static final int FIELD_AGE = FIELD_GENDER + 1;
	public static final int FIELD_ENTER_POP_AT_AGE = FIELD_AGE + 1;
	public static final int FIELD_ENTER_POP_AT_TIME = FIELD_ENTER_POP_AT_AGE + 1;
	public static final int FIELD_MAX_REGULAR_PARTNER_12_MONTHS = FIELD_ENTER_POP_AT_TIME + 1;
	public static final int FIELD_MAX_CASUAL_PARTNERS_12_MONTHS = FIELD_MAX_REGULAR_PARTNER_12_MONTHS + 1;

	// Int[] field
	public static final int FIELD_INFECT_STAT = FIELD_MAX_CASUAL_PARTNERS_12_MONTHS + 1;
	public static final int FIELD_LAST_INFECTED_AT_AGE = FIELD_INFECT_STAT + 1;
	public static final int FIELD_TIME_UNTIL_NEXT_STAGE = FIELD_LAST_INFECTED_AT_AGE + 1;
	public static final int FIELD_NUM_PARTNER = FIELD_TIME_UNTIL_NEXT_STAGE + 1;

	// boolean[] field
	public static final int FIELD_LAST_ACT_INFECTIOUS = FIELD_NUM_PARTNER + 1;

	private static final int FIELDTYPE_LIMIT_INT = FIELD_INFECT_STAT;
	private static final int FIELDTYPE_LIMIT_INT_ARRAYS = FIELD_LAST_ACT_INFECTIOUS;

	public static final int GENDER_TYPE_FEMALE = 0;
	public static final int GENDER_TYPE_HETRO_MALE = 1;
	public static final int GENDER_TYPE_MSMO = 2;
	public static final int GENDER_TYPE_MSMW = 3;

	private final Pattern PARAM_PATTERN = Pattern.compile("\\d+");

	// Casual encounter record
	protected int[] casualRecord = new int[ONE_YEAR_INT];
	protected int casualRecordIndex = 0;
	
	// Regular partner record - 12 months
	protected int[] regularRecord = new int[ONE_YEAR_INT];
	protected int regularRecordIndex = 0;

	public Person_Bridging_Pop(int id, int gender, double startingAge, int startingTime, int numInfection) {
		super();
		fields[FIELD_ID] = id;
		fields[FIELD_GENDER] = gender;
		fields[FIELD_AGE] = startingAge;
		fields[FIELD_ENTER_POP_AT_AGE] = startingAge;
		fields[FIELD_ENTER_POP_AT_TIME] = startingTime;

		for (int i = FIELDTYPE_LIMIT_INT; i < FIELDTYPE_LIMIT_INT_ARRAYS; i++) {
			fields[i] = new int[numInfection];
			Arrays.fill((int[]) fields[i], AbstractIndividualInterface.INFECT_S);
		}

		fields[FIELD_LAST_ACT_INFECTIOUS] = new boolean[numInfection];
	}

	public void addCasualPartner(AbstractIndividualInterface p) {
		casualRecord[casualRecordIndex] = p.getId();
	}
	
	public void addRegularPartner(AbstractIndividualInterface p) {
		regularRecord[regularRecordIndex] = p.getId();
	}
	
	public int getNumRegularInRecrod() {
		return getNumInRecord(regularRecord);
	}
	
	public int getNumCasualInRecord() {
		return getNumInRecord(casualRecord);
	}
	
	private int getNumInRecord(int [] record) {
		int c = 0;
		for(int i = 0; i < record.length; i++) {
			if(record[i] > 0) {
				c++;
			}
		}
		return c;
	}


	public int getGenderType() {
		return (int) fields[FIELD_GENDER];
	}

	@Override
	public double getAge() {
		return (double) fields[FIELD_AGE];
	}

	@Override
	public int getId() {
		return (int) fields[FIELD_ID];
	}

	@Override
	public int[] getInfectionStatus() {
		return (int[]) fields[FIELD_INFECT_STAT];
	}

	@Override
	public int getInfectionStatus(int index) {
		return getInfectionStatus()[index];
	}

	@Override
	public double getLastInfectedAtAge(int infectionIndex) {
		return ((int[]) fields[FIELD_LAST_INFECTED_AT_AGE])[infectionIndex];
	}

	@Override
	public Comparable<?> setParameter(String id, Comparable<?> value) {
		Matcher m = PARAM_PATTERN.matcher(id);
		Comparable<?> res = null;

		if (m.find()) {
			int nId = Integer.parseInt(m.group());

			if (nId < FIELDTYPE_LIMIT_INT) {
				res = (int) fields[nId];
				fields[nId] = (int) value;
			} else if (nId < FIELDTYPE_LIMIT_INT_ARRAYS) {
				if (m.find()) {
					int aInd = Integer.parseInt(m.group());
					res = ((int[]) fields[nId])[aInd];
					((int[]) fields[nId])[aInd] = (int) value;
				}

			} else {
				if (m.find()) {
					int aInd = Integer.parseInt(m.group());
					res = ((boolean[]) fields[nId])[aInd];
					((boolean[]) fields[nId])[aInd] = (boolean) value;
				}

			}

		}

		if (res == null) {
			System.err.println(String.format("setParameter: Ill-formed parameter string of \"%s\" and value \"%s\"", id,
					value.toString()));
		}
		return res;
	}

	@Override
	public Comparable<?> getParameter(String id) {
		Matcher m = PARAM_PATTERN.matcher(id);
		if (m.find()) {
			int nId = Integer.parseInt(m.group());

			if (nId < FIELDTYPE_LIMIT_INT) {
				return (int) fields[nId];
			} else if (nId < FIELDTYPE_LIMIT_INT_ARRAYS) {
				if (m.find()) {
					int aInd = Integer.parseInt(m.group());
					return ((int[]) fields[nId])[aInd];
				}
			} else {
				if (m.find()) {
					int aInd = Integer.parseInt(m.group());
					return ((boolean[]) fields[nId])[aInd];
				}
			}

		}

		System.err.println(String.format("getParameter: Ill-formed parameter string of \"%s\" ", id));
		return null;
	}

	@Override
	public double getTimeUntilNextStage(int index) {
		return ((int[]) fields[FIELD_TIME_UNTIL_NEXT_STAGE])[index];
	}

	@Override
	public boolean isMale() {
		return ((int) fields[FIELD_GENDER]) != 0; // Include hetro male, MSMO and MSMW
	}

	@Override
	public void setAge(double age) {
		fields[FIELD_AGE] = (int) age;
	}

	@Override
	public void setInfectionStatus(int index, int newInfectionStatus) {
		((int[]) fields[FIELD_INFECT_STAT])[index] = newInfectionStatus;
	}

	@Override
	public void setLastActInfectious(int infectionIndex, boolean lastActInf) {
		((boolean[]) fields[FIELD_LAST_ACT_INFECTIOUS])[infectionIndex] = lastActInf;

	}

	@Override
	public void setTimeUntilNextStage(int index, double newTimeUntilNextStage) {
		((int[]) fields[FIELD_TIME_UNTIL_NEXT_STAGE])[index] = (int) newTimeUntilNextStage;
	}

	@Override
	public int getEnterPopulationAt() {
		return (int) fields[FIELD_ENTER_POP_AT_TIME];
	}

	@Override
	public double getStartingAge() {
		return (int) fields[FIELD_ENTER_POP_AT_AGE];
	}

	@Override
	public void setEnterPopulationAt(int enterPopulationAt) {
		fields[FIELD_ENTER_POP_AT_TIME] = enterPopulationAt;

	}

	@Override
	public void setLastInfectedAtAge(int infectionIndex, double age) {
		((int[]) fields[FIELD_LAST_INFECTED_AT_AGE])[infectionIndex] = (int) age;

	}

	@Override
	public int incrementTime(int deltaT, AbstractInfection[] infectionList) {
		this.setAge(this.getAge() + deltaT);
	
		int res = deltaT;
		for (int t = 0; t < deltaT; t++) {			
			
			casualRecordIndex = (casualRecordIndex + 1) % casualRecord.length;
			if (casualRecord[casualRecordIndex] != 0) {
				// Remove casual from record
				casualRecord[casualRecordIndex] = 0;
			}
			
			regularRecordIndex = (regularRecordIndex + 1) % regularRecord.length;
			if (regularRecord[regularRecordIndex] != 0) {
				// Remove regular from record
				regularRecord[regularRecordIndex] = 0;
			}
		}
		return res;
	}

}
