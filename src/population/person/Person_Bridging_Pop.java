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
			18 * AbstractIndividualInterface.ONE_YEAR_INT,
			// FIELD_ENTER_POP_AT_AGE
			18 * AbstractIndividualInterface.ONE_YEAR_INT,
			// FIELD_ENTER_POP_AT
			0,
			// FIELD_INFECT_STAT
			new int[0],
			// FIELD_LAST_INFECTED_AT_AGE
			new int[0],
			// FIELD_TIME_UNTIL_NEXT_STAGE
			new int[0],
			// FIELD_LAST_ACT_INFECTIOUS
			new boolean[0],

	};

	public static final int FIELD_ID = 0;
	public static final int FIELD_GENDER = FIELD_ID + 1;
	public static final int FIELD_AGE = FIELD_GENDER + 1;
	public static final int FIELD_ENTER_POP_AT_AGE = FIELD_AGE + 1;
	public static final int FIELD_ENTER_POP_AT_TIME = FIELD_ENTER_POP_AT_AGE + 1;
	public static final int FIELD_INFECT_STAT = FIELD_ENTER_POP_AT_TIME + 1;
	public static final int FIELD_LAST_INFECTED_AT_AGE = FIELD_INFECT_STAT + 1;
	public static final int FIELD_TIME_UNTIL_NEXT_STAGE = FIELD_LAST_INFECTED_AT_AGE + 1;
	public static final int FIELD_LAST_ACT_INFECTIOUS = FIELD_TIME_UNTIL_NEXT_STAGE + 1;
	public static final int FIELD_BEHAVIOR_TYPE = FIELD_LAST_ACT_INFECTIOUS + 1;
	
	
	public static final int GENDER_TYPE_FEMALE = 0;
	public static final int GENDER_TYPE_HETRO_MALE = 1;
	public static final int GENDER_TYPE_MSMO = 2;
	public static final int GENDER_TYPE_MSMW = 3;
	
	public static final int BEHAVIOR_TYPE_REGULAR_ONLY = 0;
	public static final int BEHAVIOR_TYPE_CASUAL_ONLY = BEHAVIOR_TYPE_REGULAR_ONLY + 1;
	public static final int BEHAVIOR_TYPE_ANY = BEHAVIOR_TYPE_CASUAL_ONLY + 1;

	private final Pattern PARAM_PATTERN = Pattern.compile("\\d+");

	public Person_Bridging_Pop(int id, int gender, int startingAge, int startingTime, int numInfection) {
		super();
		fields[FIELD_ID] = id;
		fields[FIELD_GENDER] = gender;
 		fields[FIELD_AGE] = startingAge;
		fields[FIELD_ENTER_POP_AT_AGE] = startingAge;
		fields[FIELD_ENTER_POP_AT_TIME] = startingTime;
		
		for(int i = FIELD_INFECT_STAT; i < FIELD_LAST_ACT_INFECTIOUS; i++) {
			fields[i] = new int[numInfection];	
			Arrays.fill((int[]) fields[i], AbstractIndividualInterface.INFECT_S);
		}
		
		fields[FIELD_LAST_ACT_INFECTIOUS] = new boolean[numInfection];
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
	public Comparable<?> getParameter(String id) {
		Matcher m = PARAM_PATTERN.matcher(id);
		if (m.find()) {
			int nId = Integer.parseInt(m.group());

			switch (nId) {
			// Int fields
			case FIELD_ID:
			case FIELD_GENDER:
			case FIELD_AGE:
			case FIELD_ENTER_POP_AT_AGE:
			case FIELD_ENTER_POP_AT_TIME:
				return (int) fields[nId];

			// Int Array
			case FIELD_INFECT_STAT:
			case FIELD_LAST_INFECTED_AT_AGE:
			case FIELD_TIME_UNTIL_NEXT_STAGE:
				if (m.find()) {
					int aInd = Integer.parseInt(m.group());
					return ((int[]) fields[nId])[aInd];
				}

			// Boolean Array
			case FIELD_LAST_ACT_INFECTIOUS:
				if (m.find()) {
					int aInd = Integer.parseInt(m.group());
					return ((boolean[]) fields[nId])[aInd];
				}
			default:
				// Nothing
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
	public Comparable<?> setParameter(String id, Comparable<?> value) {
		Matcher m = PARAM_PATTERN.matcher(id);
		Comparable<?> res = null;

		if (m.find()) {
			int nId = Integer.parseInt(m.group());

			switch (nId) {
			// Int fields
			case FIELD_ID:
			case FIELD_GENDER:
			case FIELD_AGE:
			case FIELD_ENTER_POP_AT_AGE:
			case FIELD_ENTER_POP_AT_TIME:
				res = (int) fields[nId];
				fields[nId] = (int) value;
			
			// Int Array
			case FIELD_INFECT_STAT:
			case FIELD_LAST_INFECTED_AT_AGE:
			case FIELD_TIME_UNTIL_NEXT_STAGE:
				if (m.find()) {
					int aInd = Integer.parseInt(m.group());
					res = ((int[]) fields[nId])[aInd];
					((int[]) fields[nId])[aInd] = (int) value;
				}

			// Boolean Array
			case FIELD_LAST_ACT_INFECTIOUS:
				if (m.find()) {
					int aInd = Integer.parseInt(m.group());
					res = ((boolean[]) fields[nId])[aInd];
					((boolean[]) fields[nId])[aInd] = (boolean) value;
				}
			default:
				// Nothing
			}

		}

		if (res == null) {
			System.err.println(String.format("setParameter: Ill-formed parameter string of \"%s\" and value \"%s\"", id,
					value.toString()));
		}
		return res;
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
		// TODO Auto-generated method stub
		return 0;
	}

}
