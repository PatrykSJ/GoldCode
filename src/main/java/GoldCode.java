import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.PrintWriter;

import java.lang.reflect.Array;
import java.util.*;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.DefaultXYDataset;


public class GoldCode {
    public StringBuilder resultToFile = new StringBuilder();

    public GoldCode() {
    }

    public class Lfsr{
        // HashMap with entire output
        private final HashMap<Integer, String> receivedChains = new HashMap<>();

        public Lfsr() {

        }

        /**
         * function that run "software lfsr"
         * @param initBitChain seed of lfsr - string without any spaces or commas
         * @param enteredPolynomial list of the degrees that taps are taken of (equal value as expected degrees in polynomial)
         * @return sequence of last bits after each "iteration/shift"
         */
        public List<Character> doLfsr( List<Integer> enteredPolynomial, String initBitChain){
            char[] currentLfsr = initBitChain.toCharArray(); // current lsfr array
            int i=1;
            double pow = Math.pow(2, initBitChain.length())-2;
            List<Character> outputLastChar = new ArrayList<>(); // list contains last character after each iteration of shift (code)

            String entireChainOutput="";
            while(true){
                char[] previousLfsr = Arrays.copyOf(currentLfsr, currentLfsr.length);
                outputLastChar.add(previousLfsr[previousLfsr.length-1]);
                byte toXor; //result with all xor operations
                toXor =  (byte) (currentLfsr[enteredPolynomial.get(0)-1] ^ currentLfsr[enteredPolynomial.get(1)-1]);
                for(int j=2; j<enteredPolynomial.size(); j++){  //xoring specified values
                    toXor =   (byte) (toXor ^  currentLfsr[enteredPolynomial.get(j)-1]);
                }
                for(int y = 1; y< currentLfsr.length; y++){ //shifting values in lsfr
                    currentLfsr[y]=previousLfsr[y-1];
                }
                if(Character.forDigit(toXor, 10)=='1'||Character.forDigit(toXor, 10)=='0'){
                    currentLfsr[0] = Character.forDigit(toXor, 10);
                }else{
                    try {
                        currentLfsr[0] = new String(new byte[]{toXor}, "ASCII").charAt(0);
                    }catch (Exception e){
                        System.out.println(e);
                    }
                }
                 //adding xor product as '0' lfsr's state
                entireChainOutput = new String(currentLfsr); //calculated current chain

                 //if the cycle repeats, break the loop
                receivedChains.put(i, entireChainOutput );
                if(i==pow) {
                    outputLastChar.add(currentLfsr[currentLfsr.length-1]);
                    break;
                }
                i++;
                if(entireChainOutput.equals(initBitChain)) {
                    System.out.println("It's not a MLS");
                    break;
                }
            }
            receivedChains.put(0, initBitChain);
            System.out.println("LFSR states: " + receivedChains);
            resultToFile.append("LFSR states: ").append(receivedChains).append("\n");
            return outputLastChar;
        }
    }

    public ArrayList<ArrayList<ArrayList<Integer>>> readMlsPolynomialsFromFile(String fileName){
        ArrayList<ArrayList<ArrayList<Integer>>> readPolynomials = new ArrayList<>();
        ArrayList<ArrayList<Integer>> fillFirstOnesWith0 = new ArrayList<>();
        ArrayList<Integer> justZero = new ArrayList<>();
        justZero.add(0);
        fillFirstOnesWith0.add(justZero);
        readPolynomials.add(fillFirstOnesWith0);
        readPolynomials.add(fillFirstOnesWith0);
        readPolynomials.add(fillFirstOnesWith0);
        readPolynomials.add(fillFirstOnesWith0);
        File myObj = new File(fileName);
        ArrayList<ArrayList<Integer>> currentPolynomials = new ArrayList<>();
        try (Scanner reader = new Scanner(myObj)){
            while (reader.hasNextLine()) {
                ArrayList<Integer> currentPolynomial = new ArrayList<>();
                String data = reader.nextLine();
                if(data.equals("//")) {
                    readPolynomials.add(currentPolynomials);
                    currentPolynomials = new ArrayList<>();
                    continue;
                }
                if(data.equals("")) break;
                int currentlyReadPolynomial = Integer.parseInt(data, 16);
                String binaryPoly = Integer.toBinaryString(currentlyReadPolynomial);
                binaryPoly = new StringBuilder(binaryPoly).reverse().toString();
                String[] ar = binaryPoly.split("");
                for (int i=0; i<ar.length; i++) {
                    if (ar[i].equals("1")) currentPolynomial.add(i+1);
                    }
                currentPolynomials.add(currentPolynomial);
            }
        } catch (Exception e) {
            System.err.println(e);
        }
        readPolynomials.add(currentPolynomials);
        return readPolynomials;
    }


    public int[] convertToArray(String s){
        int[] t = new int[s.length()];
        int j=0;
        for(int i=s.length()-1; i>=0; i--) {
            t[j] = Integer.parseInt(String.valueOf(s.charAt(i))); //od lewej najmniejsze potegi
            j++;
        }
        return t;
    }


    public boolean algorithmCheck(ArrayList<ArrayList<ArrayList<Integer>>> allPolynomialsMls, ArrayList<Integer> firstPolynomialCoff, ArrayList<Integer> secondPolynomialCoff, String polynomial){ // najmniejszy element na poczatku list
        int n = polynomial.length()-1;
        double l;
        //if(firstPolynomialCoff.get(0)%2!=0){
        if(n%2!=0){
            l = Math.pow(2, firstPolynomialCoff.get(0))+1;
        }else{
            l = Math.pow(2, (firstPolynomialCoff.get(0)+2)/2.)+1;
        }
        ArrayList<ArrayList<Integer>> nPolynomials = allPolynomialsMls.get(n);
        for(int i=0 ;i<nPolynomials.size(); i++){
            if(secondPolynomialCoff.get(0) != l) return false;
            if(nPolynomials.get(i).get(0) == l && nPolynomials.get(i).equals(secondPolynomialCoff)){
                double level;
                if(n%2==0){
                    level = Math.pow(2, (n+1)/2.) +1;
                }else{
                    level = Math.pow(2, (n+1)/2.) -1;
                }

                double db;
                db = 20*Math.log10((Math.pow(2, n))/level);
                System.out.println("Correlation level: " + level);
                System.out.println("Difference in [dB]: " + db);
                return true;
            }
        }
        return false;
    }

    public boolean checkAllConditions(String polynomial1, String polynomial2, ArrayList<Integer> firstPolynomialCoff, ArrayList<Integer> secondPolynomialCoff, ArrayList<ArrayList<ArrayList<Integer>>> allPolynomialsMls){
        if(polynomial1.length()!=polynomial2.length()) return false;
        return algorithmCheck(allPolynomialsMls, firstPolynomialCoff, secondPolynomialCoff, polynomial1);
    }
    /**
     * Function GC calculates gold code's value
     * @param enteredPolynomial first polynomial - list of integers (degrees which taps values are taken to xor)
     * @param enteredPolynomial2 second polynomial - same rules applied
     * @param initBitChain string of initial seed for the second polynomial (without any spaces or commas)
     * @return gold code as List of Integers
     */
    private List<Integer> GC(List<Integer> enteredPolynomial, List<Integer> enteredPolynomial2, String initBitChain) {
        List<Integer> goldCode = new ArrayList<>();

        Lfsr lfsr = new Lfsr();

        List<Character> goldcode1, goldcode2 = new ArrayList<>();
        StringBuilder defInitChain = new StringBuilder(); //building default seed for first polynomial
        for(int i =0 ; i<initBitChain.length()-1; i++){
            defInitChain.append(0);
        }
        defInitChain.append(1);
        goldcode1 = lfsr.doLfsr(enteredPolynomial, defInitChain.toString()); //calculating lfsr's output for the first polynomial
        System.out.println("Output of the first lfsr:" + goldcode1);
        resultToFile.append("Output of the first lfsr:").append(goldcode1).append("\n");

        goldcode2 = lfsr.doLfsr(enteredPolynomial2, initBitChain);
        System.out.println("Output of the second lfsr:" + goldcode2); //calculating lfsr's output for the second polynomial
        resultToFile.append("Output of the second lfsr:").append(goldcode2).append("\n");

        for(int i=0; i< goldcode1.size(); i++){
            goldCode.add(goldcode1.get(i) ^ goldcode2.get(i));  //xor of the two lfsrs' output
        }


        System.out.println("Final Gold Code value: " + goldCode); //printing and returning final gold code value
        resultToFile.append("Final Gold Code value: ").append(goldCode).append("\n");
        return goldCode;
    }

    private List<Integer> GCwithItself(List<Integer> enteredPolynomial, String initBitChain) {
        List<Integer> goldCode = new ArrayList<>();
        Lfsr lfsr = new Lfsr();

        List<Character> goldcode1 = new ArrayList<>();
        StringBuilder defInitChain = new StringBuilder(); //building default seed for first polynomial
        for(int i =0 ; i<initBitChain.length()-1; i++){
            defInitChain.append(0);
        }
        defInitChain.append(1);
        goldcode1 = lfsr.doLfsr(enteredPolynomial, defInitChain.toString()); //calculating lfsr's output for the first polynomial
        System.out.println("Output of the first lfsr:" + goldcode1);
        resultToFile.append("Output of the first lfsr:").append(goldcode1).append("\n");

        for(int i=0; i< goldcode1.size(); i++){
            goldCode.add((int)goldcode1.get(i));  //xor of the two lfsrs' output
        }

        System.out.println("Final Gold Code value: " + goldcode1); //printing and returning final gold code value
        resultToFile.append("Final Gold Code value: ").append(goldcode1).append("\n");
        return goldCode;
    }

    /**
     * Funtion calculates autocorrelation
     * @param givenCode given code as list of integers
     */
    private void calculateAutoCorrelation(List<Integer> givenCode){
        List<Integer> calculatedACor = new ArrayList<>();
        List<Integer> shiftedCode = new ArrayList<>(givenCode);
        for (int i=0; i< givenCode.size(); i++){
            int sum =0;                                 //d distance value
            for(int j=0; j< givenCode.size(); j++){
                if(givenCode.get(j).equals(shiftedCode.get(j))) sum +=1; //checking if corresponding values are equal
                else sum -= 1;
            }
            calculatedACor.add(sum);
            List<Integer> previousShiftedCode = new ArrayList<>(shiftedCode);
            Integer copyOfLastEl = shiftedCode.get(shiftedCode.size()-1);

            for(int y = 1; y< previousShiftedCode.size(); y++){
                shiftedCode.set(y, previousShiftedCode.get(y-1));  //shifting chain of bits
            }
            shiftedCode.set(0, copyOfLastEl);
        }
        System.out.println("Calculated autocorrelation: " + calculatedACor); // printing calculated values
        drawChart(calculatedACor, "Unnormalized autocorrelation");
        resultToFile.append("Calculated autocorrelation: ").append(calculatedACor).append("\n");
    }


    private void calculateNormalizedAutoCorrelation(List<Integer> givenCode, String polynomial){
        List<Integer> calculatedACor = new ArrayList<>();
        List<Integer> shiftedCode = new ArrayList<>(givenCode);
        Collections.replaceAll(calculatedACor,0,-1);
        Collections.replaceAll(shiftedCode,0,-1);
        for (int i=0; i< givenCode.size(); i++){
            int sum =0;                                 //d distance value
            for(int j=0; j< givenCode.size(); j++){
                sum += givenCode.get(j)*shiftedCode.get(j);
            }
            if(i==0){
                calculatedACor.add(sum / (int) Math.pow(2, (polynomial.length() - 1))-1);
            }else {
                calculatedACor.add(sum / (polynomial.length() - 1));
            }
            List<Integer> previousShiftedCode = new ArrayList<>(shiftedCode);
            Integer copyOfLastEl = shiftedCode.get(shiftedCode.size()-1);

            for(int y = 1; y< previousShiftedCode.size(); y++){
                shiftedCode.set(y, previousShiftedCode.get(y-1));  //shifting chain of bits
            }
            shiftedCode.set(0, copyOfLastEl);
        }
        System.out.println("Calculated autocorrelation: " + calculatedACor); // printing calculated values
        drawChart(calculatedACor, "Normalized autocorrelation");
        resultToFile.append("Calculated autocorrelation: ").append(calculatedACor).append("\n");
    }


    public static String convertToString(StringBuilder sb){
        return sb.toString();
    }

    public void drawChart(List<Integer> correlation, String title){
        double[] from0ToP = new double[correlation.size()];
        double[]  valuesCorrelation =  new double[correlation.size()];
        for(int i=0; i<correlation.size(); i++){
            from0ToP[i] = (double) i;
            valuesCorrelation[i] = correlation.get(i);
        }
        DefaultXYDataset dataset = new DefaultXYDataset();
        double[][] data = new double[][] {from0ToP, valuesCorrelation};
        dataset.addSeries("AutoCorrelation", data);
        JFreeChart chart = ChartFactory.createScatterPlot(title, "Sequence Index", "AutoCorrelation Values", dataset);
        ChartFrame frame = new ChartFrame("Chart", chart);
        frame.pack();
        frame.setVisible(true);
        ChartPanel panel = new ChartPanel(chart);
        panel.setPreferredSize(new java.awt.Dimension(500, 300));
    }

    public ArrayList<Integer> calculatePossibleValues(String polynomial){
        ArrayList<Integer> values = new ArrayList<>();
        int n = polynomial.length()-1;
        if(n%2!=0){
            int t = (int) (1+Math.pow(2, (n+1)/2.));
            values.add(-t);
            values.add(-1);
            values.add(t-2);
        }else{
            int t = (int) (1+Math.pow(2, (n+2)/2.));
            values.add(-t);
            values.add(-1);
            values.add(t-2);
        }
        return values;
    }


    public void writeSeedsToFile(String file, int n){
        try ( PrintWriter out = new PrintWriter(file)){
            int pow = (int) Math.pow(2, n);
            for(int i=1; i<pow; i++){
                StringBuilder sb = new StringBuilder();
                String afterconvert = Integer.toBinaryString(i);
                if(afterconvert.length()<n){
                    for(int j=0; j<n-afterconvert.length(); j++){
                        sb.append("0");
                    }
                }
                sb.append(afterconvert);
                out.println(sb.toString());
            }
        }catch(Exception e){
            System.out.println(e);
        }
    }


    public void calculateCrossCorelation(List<Integer> enteredPolynomial, List<Integer> enteredPolynomial2, int n, String file){
        writeSeedsToFile(file, n);
        List<Integer> goldCodeConst = new ArrayList<>();
        List<Integer> goldCodePol2 = new ArrayList<>();
        Lfsr lfsr = new Lfsr();

        List<Character> goldcode1, goldcode2, gold2Const = new ArrayList<>(); //gold2Cobst nie potrzebne
        StringBuilder defInitChain = new StringBuilder(); //building default seed for first polynomial
        for(int i =0 ; i<n-1; i++){
            defInitChain.append(0);
        }
        defInitChain.append(1); //const
        goldcode1 = lfsr.doLfsr(enteredPolynomial, defInitChain.toString()); //calculating lfsr's output for the first polynomial

        //goldcode for 2nd generation
        String input;
        for (int i = 0; i < goldcode1.size(); i++) {
            goldCodeConst.add(goldcode1.get(i) - '0');
        }
        // do wywalenia w razie cross LFSR
        ArrayList<Integer> defaultGoldCode = new ArrayList<>();
        gold2Const = lfsr.doLfsr(enteredPolynomial2, defInitChain.toString());
        for(int i=0; i< goldcode1.size(); i++){
            defaultGoldCode.add(goldcode1.get(i) ^ gold2Const.get(i));  //xor of the two lfsrs' output
        }

        File myObj = new File(file);
        try(Scanner fileReader = new Scanner(myObj)){
            while(fileReader.hasNextLine()) {
                input = fileReader.nextLine();
                goldcode2 = lfsr.doLfsr(enteredPolynomial2, input);
                for (int i = 0; i < goldcode2.size(); i++) {
                    goldCodePol2.add(goldcode2.get(i) - '0');
                }
                // do wywalenia
                ArrayList<Integer> TempGoldCode = new ArrayList<>();
                for(int i=0; i< goldcode1.size(); i++){
                    TempGoldCode.add(goldcode1.get(i) ^ goldcode2.get(i));  //xor of the two lfsrs' output
                }

                double exp;
                if(n%2==0){
                    exp = Math.pow(2, (n+2)/2.)+1;
                }else{
                    exp = Math.pow(2, (n+1)/2.)+1;
                }


                // Cross correlacja kodow golda
                List<Integer> calculatedACor = new ArrayList<>();
                for (int i = 0; i < defaultGoldCode.size(); i++) {
                    int sum = 0;
                    for (int j = 0; j < defaultGoldCode.size(); j++) {
                        if (defaultGoldCode.get(j).equals(TempGoldCode.get(j))) sum += 1;
                        else sum -= 1;
                    }
                    calculatedACor.add(sum);
                    double calc = Math.pow(2, n)-1;
                    if(Math.abs(sum) > exp&&Math.abs(sum)!=calc){
                        System.out.println(sum);
                        System.out.println("incorrect autocorrelation");
                        return;
                    }
                    List<Integer> previousShiftedCode = new ArrayList<>(TempGoldCode);
                    Integer copyOfLastEl = TempGoldCode.get(TempGoldCode.size() - 1);

                    for (int y = 1; y < previousShiftedCode.size(); y++) {
                        TempGoldCode.set(y, previousShiftedCode.get(y - 1));
                    }
                    TempGoldCode.set(0, copyOfLastEl);
                }

                System.out.println(calculatedACor);

               /*
                //cross correlation
                List<Integer> calculatedACor = new ArrayList<>();
                for (int i = 0; i < goldCodeConst.size(); i++) {
                    int sum = 0;
                    for (int j = 0; j < goldCodeConst.size(); j++) {
                        if (goldCodeConst.get(j).equals(goldCodePol2.get(j))) sum += 1;
                        else sum -= 1;
                    }
                    calculatedACor.add(sum);
                    List<Integer> previousShiftedCode = new ArrayList<>(goldCodePol2);
                    Integer copyOfLastEl = goldCodePol2.get(goldCodePol2.size() - 1);

                    for (int y = 1; y < previousShiftedCode.size(); y++) {
                        goldCodePol2.set(y, previousShiftedCode.get(y - 1));
                    }
                    goldCodePol2.set(0, copyOfLastEl);
                }
                System.out.println(calculatedACor);

*/

            }
        }catch (Exception e){
            System.out.println(e);
        }
    }


    public static void main(String[] args) {
        //https://www.gaussianwaves.com/2015/06/gold-code-generator/
        //https://www.mathworks.com/help/comm/ref/comm.goldsequence-system-object.html#mw_e154c06f-3458-4b1b-b6de-73dca95545a1
        ArrayList<ArrayList<ArrayList<Integer>>> readPolynomialsMls = new ArrayList<>();
        GoldCode GcClient = new GoldCode();
        readPolynomialsMls = GcClient.readMlsPolynomialsFromFile("mls_back.txt");
        File myObj = new File("input.txt");
        int line=0;
            ArrayList<Integer> list1 = new ArrayList<>();
            ArrayList<Integer> list2 = new ArrayList<>();
            String pol1="", pol2="";
            try {
                String seedForSecondPolynomial ="";
                    try(Scanner fileReader = new Scanner(myObj)){
                        while (fileReader.hasNextLine()) {
                            String data, input;
                            input = fileReader.nextLine();
                            if(line%3==0) pol1 = input;
                            else if (line%3==1) pol2 = input;
                            if(line%3==2) data = input; //od lewej najwyzsza potega
                            else data = new StringBuilder(input).reverse().toString();
                            String[] ar = data.split("");
                            for (int i=1; i<ar.length; i++) {
                                if (line % 3 == 0 && ar[i].equals("1"))list1.add(i);
                                else if (line % 3 == 1 && ar[i].equals("1")) list2.add(i);
                                else if(line%3==2) {
                                    seedForSecondPolynomial = data;
                                    break;
                                }
                            }
                            line++;

                        }
                    }catch (Exception e){
                        System.out.println(e);
                    }
                    if(list2.isEmpty()){
                        Collections.reverse(list1);
                        List<Integer> calculatedGoldCode = GcClient.GCwithItself(list1, seedForSecondPolynomial);
                        GcClient.calculateAutoCorrelation(calculatedGoldCode);
                    }else {
                        System.out.println(GcClient.checkAllConditions(pol1, pol2, list1, list2, readPolynomialsMls)); //tutaj listy nie zreversowane, wiec najmniejszy element na poczatku
                        System.out.println("Possible values in AutoCorrelation:" + GcClient.calculatePossibleValues(pol1));
                        Collections.reverse(list1);
                        Collections.reverse(list2);
                        List<Integer> calculatedGoldCode = GcClient.GC(list1, list2, seedForSecondPolynomial);
                        GcClient.calculateAutoCorrelation(calculatedGoldCode);
                        GcClient.calculateNormalizedAutoCorrelation(calculatedGoldCode, pol1);
                        System.out.println("\n");
                        GcClient.calculateCrossCorelation(list1, list2, 9, "seeds9N.txt");
                    }
                    try ( PrintWriter out = new PrintWriter("fileWithResults.txt")){
                        out.print(convertToString(GcClient.resultToFile));
                    }catch(Exception e){
                        System.out.println(e);
                    }
                }
            catch (Exception e){
                System.err.println(e);
            }
    }
}
