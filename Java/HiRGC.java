import java.io.*;
import java.util.*;
import java.util.Date;
import java.lang.management.*;

public class Algorithm {

    static final byte[] code_rule = {0, 1, 2, 3};
    static final char[] invert_code_rule = {'A', 'C', 'G', 'T'};
    static final int SUB_STR_LEN = 20;
    static final int MAX_CHAR_NUM = 1<<28;//1<<30
    static final int MAX_ARR_NUM_BIT = 30;
    static final int MAX_ARR_NUM_BIT_SHIFT = MAX_ARR_NUM_BIT>>1;
    static final int MAX_ARR_NUM = 1<<MAX_ARR_NUM_BIT;
    public static byte[] ref_seq_code = new byte[MAX_CHAR_NUM];
    public static byte[] tar_seq_code = new byte[MAX_CHAR_NUM];
    public static String metadata;
    public static int ref_seq_len = 0, tar_seq_len = 0;
    public static int pos_vec_len = 0, line_break_len = 0, other_char_len = 0, n_vec_len = 0;

    // public static POSITION_RANGE[] pos_vec = new POSITION_RANGE[1<<20];
    // public static POSITION_RANGE[] n_vec = new POSITION_RANGE[1<<20];
    // public static POSITION_OTHER_CHAR[] other_char_vec = new POSITION_OTHER_CHAR[1<<20];
    public static int[] line_break_vec = new int[1<<25];

    public static int[] pos_vec_begin = new int[1<<20];
    public static int[] pos_vec_length = new int[1<<20];
    public static int[] n_vec_begin = new int[1<<20];
    public static int[] n_vec_length = new int[1<<20];
    public static int[] other_char_vec_pos = new int[1<<20];
    public static byte[] other_char_vec_ch = new byte[1<<20];

    // public static int[] loc_v = new int[MAX_CHAR_NUM];
    public static int[] loc_next = new int[MAX_CHAR_NUM];
    public static int[] point = new int[MAX_ARR_NUM];

    public Algorithm() throws IOException {
        // result_text = "";
    }

    public static byte agctIndex(char ch) {
        if (ch == 'A') {
            return 0;
        }
        if (ch == 'C') {
            return 1;
        } 
        if (ch == 'G') {
            return 2;
        }
        if (ch == 'T') {
            return 3;
        }
        return -1;
    }

    public static void readRefFile(String refFile) throws IOException {
        FileInputStream fileInputStream = new FileInputStream(refFile);
        DataInputStream dataInputStream = new DataInputStream(fileInputStream);
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(dataInputStream));
        String line;
        char temp_ch;
        byte index = -1;
        int line_length;

        ref_seq_len = 0;
        bufferedReader.readLine();//first line is identifier
        while ((line = bufferedReader.readLine()) != null) {
            // System.out.printf("%s\n", line);
            line_length = line.length();
            for (int i = 0; i < line_length; i++) {
                temp_ch = line.charAt(i);
                if (!Character.isUpperCase(temp_ch)) {
                    temp_ch = Character.toUpperCase(temp_ch);
                }
                index = agctIndex(temp_ch);
                if (index > -1) {
                    ref_seq_code[ref_seq_len++] = code_rule[index];
                }
            } 
        }
        // System.out.printf("ref_seq_len: %d\n", ref_seq_len);
        // for (int i = 0; i < ref_seq_len; i++) {
        //     System.out.printf("%d", ref_seq_code[i]);
        // }
    }

    public static void readTarFile(String tarFile) throws IOException {
        FileInputStream fileInputStream = new FileInputStream(tarFile);
        DataInputStream dataInputStream = new DataInputStream(fileInputStream);
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(dataInputStream));
        String line;
        char temp_ch;
        byte index = -1;
        int line_length;
        boolean flag = true, n_flag = false;
        int letters_len = 0, n_letters_len = 0;

        pos_vec_len = 0;
        line_break_len = 0;
        other_char_len = 0;
        n_vec_len = 0;
        tar_seq_len = 0;

        metadata = bufferedReader.readLine();
        while ((line = bufferedReader.readLine()) != null) {
            // System.out.printf("%s\n", line);
            line_length = line.length();
            for (int i = 0; i < line_length; i++) {
                temp_ch = line.charAt(i);
                if (Character.isLowerCase(temp_ch)) {
                    if (flag) {
                        flag = false;
                        pos_vec_begin[pos_vec_len] = letters_len;
                        letters_len = 0;
                    }
                    temp_ch = Character.toUpperCase(temp_ch);
                } else {
                    if (Character.isUpperCase(temp_ch)) {
                        if (!flag) {
                            flag = true;
                            pos_vec_length[pos_vec_len] = letters_len;
                            pos_vec_len++;
                            letters_len = 0;
                        }       
                    }
                }
                letters_len++;

                if (temp_ch != 'N') {
                    index = agctIndex(temp_ch);
                    if (index > -1) {
                        tar_seq_code[tar_seq_len++] = code_rule[index];
                    } else {
                        other_char_vec_pos[other_char_len] = tar_seq_len;
                        other_char_vec_ch[other_char_len] = (byte)(temp_ch - 'A');
                        other_char_len++;
                    }
                }

                if (!n_flag) {
                    if (temp_ch == 'N') {
                        n_vec_begin[n_vec_len] = n_letters_len;
                        n_letters_len = 0;
                        n_flag = true;
                    }
                } else {
                    if (temp_ch != 'N') {
                        n_vec_length[n_vec_len] = n_letters_len;
                        n_vec_len++;
                        n_letters_len = 0;
                        n_flag = false;
                    }
                }

                n_letters_len++;
            }
            line_break_vec[line_break_len++] = line_length;
        }
        if (!flag) {
            pos_vec_length[pos_vec_len] = letters_len;
            pos_vec_len++;
        }

        if (n_flag) {
            n_vec_length[n_vec_len] = n_letters_len;
            n_vec_len++;
        }

        for (int i = other_char_len-1; i > 0; i--) {
            other_char_vec_pos[i] -= other_char_vec_pos[i-1];
        }
        // System.out.printf("tar_seq_len: %d\n", tar_seq_len);
        // for (int i = 0; i < tar_seq_len; i++) {
        //     System.out.printf("%d", tar_seq_code[i]);
        // }
    }

    public static long getCpuTime() {
        ThreadMXBean bean = ManagementFactory.getThreadMXBean();
        return bean.isCurrentThreadCpuTimeSupported() ? bean.getCurrentThreadCpuTime() : 0L;
    }

    public static void write(String fileName, String text, boolean append) throws IOException {
        FileWriter fileWriter = new FileWriter(fileName, append);
        BufferedWriter out = new BufferedWriter(fileWriter);
        out.write(text);
        out.close();
    }

    public static void saveOtherData(String fileName) throws IOException {
        StringBuilder strBuf = new StringBuilder(1<<25);
        strBuf.append(metadata+'\n');

        //write line_break information (always line_break_len > 0)
        StringBuilder lineBuf = new StringBuilder(1<<15);
        lineBuf.append(line_break_vec[0]);
        int pre_value = line_break_vec[0], cnt = 1, code_len = 1;
        for (int i = 1; i < line_break_len; i++) {
            if (line_break_vec[i] == pre_value) {
                cnt++;
            } else {
                lineBuf.append(' '); lineBuf.append(cnt);
                lineBuf.append(' '); lineBuf.append(line_break_vec[i]);
                pre_value = line_break_vec[i];
                code_len += 2;
                cnt = 1;
            }
        }
        lineBuf.append(' '); lineBuf.append(cnt); lineBuf.append('\n');
        code_len++;
        strBuf.append(code_len); strBuf.append(' ');
        strBuf.append(lineBuf);
        //-----------

        //write pos_vec information
        strBuf.append(pos_vec_len);
        for (int i = 0; i < pos_vec_len; i++) {
            strBuf.append(' '); strBuf.append(pos_vec_begin[i]);
            strBuf.append(' '); strBuf.append(pos_vec_length[i]);
        }
        strBuf.append('\n');
        //write n_pos_vec information
        strBuf.append(n_vec_len);
        for (int i = 0; i < n_vec_len; i++) {
            strBuf.append(' '); strBuf.append(n_vec_begin[i]); 
            strBuf.append(' '); strBuf.append(n_vec_length[i]);
        }
        strBuf.append('\n');
        //write other_char information
        strBuf.append(other_char_len);
        if (other_char_len > 0) {
            int[] flag = new int[30];
            int[] arr = new int[100];
            int arr_len = 0;

            for (int i = 0; i < 26; i++) {
                flag[i] = -1;
            }

            for (int i = 0; i < other_char_len; i++) {
                strBuf.append(' '); strBuf.append(other_char_vec_pos[i]);
                byte temp = other_char_vec_ch[i];
                if (flag[temp] == -1) {
                    arr[arr_len] = temp;
                    flag[temp] = arr_len;
                    arr_len++;
                }
            }

            //save other char
            strBuf.append(' '); strBuf.append(arr_len);
            for (int i = 0; i < arr_len; i++) {
                strBuf.append(' '); strBuf.append(arr[i]);
            }
            if (arr_len < 10) {
                strBuf.append(' ');
                for (int i = 0; i < other_char_len; i++) {
                    strBuf.append(flag[other_char_vec_ch[i]]);
                }
            } else {
                //because all < 10 xxx
            }
        }
        strBuf.append('\n');
        write(fileName, strBuf.toString(), false);
    }

    public static void preProcessRef() {
        for (int i = 0; i < MAX_ARR_NUM; i++) {
            point[i] = -1;
        }
        long value = 0;
        for (int k = SUB_STR_LEN - 1; k >= 0; k--) {
            value <<= 2;
            value += ref_seq_code[k];
        }
        int id = (int)(value&(long)(MAX_ARR_NUM - 1));
        loc_next[0] = point[id];
        point[id] = 0;

        int step_len = ref_seq_len - SUB_STR_LEN + 1;
        int shift_bit_num = SUB_STR_LEN*2-2;
        int one_sub_str = SUB_STR_LEN - 1;

        for (int i = 1; i < step_len; i++) {
            value >>= 2;
            value += ((long)ref_seq_code[i + one_sub_str]<<shift_bit_num);

            id = (int)(value&(long)(MAX_ARR_NUM - 1));
            loc_next[i] = point[id];
            point[id] = i;
        }
    }

    public static void searchMatch(String resFile) throws IOException  {
        preProcessRef();

        StringBuilder strBuf = new StringBuilder(1<<30);
        int pre_dismatch_point = 0, pre_pos = 0;
        int step_len = tar_seq_len - SUB_STR_LEN + 1;

        for (int i =0; i < step_len; i++) {
            long tar_value = 0;
            for (int k = SUB_STR_LEN - 1; k >= 0; k--) {
                tar_value <<= 2;
                tar_value += tar_seq_code[i+k];
            }

            int id = point[(int)(tar_value&(long)(MAX_ARR_NUM - 1))];
            if (id > -1) {
                int tar_int_value = (int)(tar_value>>MAX_ARR_NUM_BIT);

                int max_length = -1, max_k = -1;
                for (int k = id; k != -1; k = loc_next[k]) {
                    int ref_idx = k + MAX_ARR_NUM_BIT_SHIFT;//loc[k].pos == k
                    int tar_idx = i + MAX_ARR_NUM_BIT_SHIFT;
                    int length = MAX_ARR_NUM_BIT_SHIFT;
                    while (ref_idx < ref_seq_len && tar_idx < tar_seq_len && ref_seq_code[ref_idx] == tar_seq_code[tar_idx]) {
                        ref_idx++;
                        tar_idx++;
                        length++;
                    }
                    if (length >= SUB_STR_LEN && length > max_length) {
                        max_length = length;
                        max_k = k;
                    }
                }

                if (max_k > -1) {
                    //first print mismatch substring
                    for(int k = pre_dismatch_point; k < i; k++) {
                        strBuf.append(tar_seq_code[k]);
                    }
                    strBuf.append('\n');

                    //then print match substring
                    int cur_pos;
                    cur_pos = max_k - pre_pos;
                    pre_pos = max_k + max_length;

                    strBuf.append(cur_pos); strBuf.append(' ');
                    strBuf.append(max_length - SUB_STR_LEN); strBuf.append('\n');

                    i += max_length;
                    pre_dismatch_point = i;
                }
            }
        }
        if (tar_seq_len > pre_dismatch_point) {
            // string dismatched_str;
            for(int k = pre_dismatch_point; k < tar_seq_len; k++) {
                strBuf.append(tar_seq_code[k]);
            }
            strBuf.append('\n');
        }

        write(resFile, strBuf.toString(), true);
    }

    public static void compressFile(String refFile, String tarFile, String resFile)  throws IOException {
        // System.out.printf("compressing %s\n", tarFile);
        readRefFile(refFile);
        readTarFile(tarFile);
        saveOtherData(resFile);
        searchMatch(resFile);
    }

/*    public static void main(String[] args) throws IOException, InterruptedException {
        // System.out.printf("%d; %c\n", code_rule[2], invert_code_rule[1]);
        // readRefFile("simu-ko-131.fa");
        // readTarFile("simu-ko-224.fa");
        // saveOtherData("ko.cd");
        // compressFile("simu-ko-131.fa", "simu-ko-224.fa", "ko.cd");
        // compressFile("../HG19/chr1.fa", "../YH/chr1.fa", "ko.cd");
        if (args.length != 3) {
            System.out.println("PLEASE, PROVIDE RIGHT NUMBER OF ARGUMENTS!");
            System.exit(0);
        }
        //run example: java -Xmx20g Algorithm ../HG18 ../KOREF20090224 KOREF20090224_ref_HG18
        String reference_fold = args[0];
        String target_fold = args[1];
        String final_fold = args[2];

        Runtime.getRuntime().exec("mkdir " + final_fold);

        Date startDate = new Date();
        long startCpuTimeNano = getCpuTime();
        System.out.println("STARTING AT: " + startDate);

        String fileName[] = new String[]{"chr1.fa", "chr2.fa", "chr3.fa", "chr4.fa", 
        "chr5.fa", "chr6.fa", "chr7.fa", "chr8.fa", "chr9.fa", "chr10.fa", 
        "chr11.fa", "chr12.fa", "chr13.fa", "chr14.fa", "chr15.fa", "chr16.fa", "chr17.fa", 
        "chr18.fa", "chr19.fa", "chr20.fa", "chr21.fa", "chr22.fa", "chrX.fa", "chrY.fa", };
        
        for (int i = 0; i < 24; i++) {
            compressFile(reference_fold+'/'+fileName[i], target_fold+'/'+fileName[i], final_fold+'/'+fileName[i]);
        }

        //./7za a YH_refHG18.7z YH_refHG18 -m0=PPMD

        String command = "./7za a " + final_fold + ".7z " + final_fold + " -m0=PPMD";
        System.out.println(command);

        String line;
        Process p = Runtime.getRuntime().exec(command);

        BufferedReader bri = new BufferedReader
                (new InputStreamReader(p.getInputStream()));

        BufferedReader bre = new BufferedReader
                (new InputStreamReader(p.getErrorStream()));

        while ((line = bri.readLine()) != null) {
            System.out.println(line);
        }
        bri.close();

        while ((line = bre.readLine()) != null) {
            System.out.println(line);
        }
        bre.close();

        p.waitFor();

        Date endDate = new Date();
        long taskCpuTimeNano = getCpuTime() - startCpuTimeNano;

        double diff = endDate.getTime() - startDate.getTime();
        double diffSeconds = diff / 1000;
        double diffMinutes = diff / (60 * 1000);

        System.out.println("Time in seconds*: " + (double) taskCpuTimeNano / 1000000000.0 + " seconds.");
        System.out.println("Time in seconds: " + diffSeconds + " seconds.");
        System.out.println("Time in minutes: " + diffMinutes + " minutes.");
    }
*/
    public static void main(String[] args) throws IOException, InterruptedException {
        //run example: java -Xmx20g Algorithm ../HG18 ../KOREF20090224 KOREF20090224_ref_HG18
        String fold_list[] = new String[]{"YH", "HG18", "HG19", "HG38", "KOREF20090131", "KOREF20090224"};

        for (int ri = 0; ri < 6; ri++) {
            for (int ti = 0; ti < 6; ti++) {
                if (ri != ti) {
                    String reference_fold = fold_list[ri];
                    String target_fold = fold_list[ti];
                    String final_fold = target_fold + "_ref_" + reference_fold;

                    Runtime.getRuntime().exec("mkdir " + final_fold);

                    Date startDate = new Date();
                    long startCpuTimeNano = getCpuTime();
                    System.out.println("STARTING AT: " + startDate);

                    String fileName[] = new String[]{"chr1.fa", "chr2.fa", "chr3.fa", "chr4.fa", 
                    "chr5.fa", "chr6.fa", "chr7.fa", "chr8.fa", "chr9.fa", "chr10.fa", 
                    "chr11.fa", "chr12.fa", "chr13.fa", "chr14.fa", "chr15.fa", "chr16.fa", "chr17.fa", 
                    "chr18.fa", "chr19.fa", "chr20.fa", "chr21.fa", "chr22.fa", "chrX.fa", "chrY.fa", };
                    
                    for (int i = 0; i < 24; i++) {
                        compressFile(reference_fold+'/'+fileName[i], target_fold+'/'+fileName[i], final_fold+'/'+fileName[i]);
                    }

                    //./7za a YH_refHG18.7z YH_refHG18 -m0=PPMD

                    String command = "./7za a " + final_fold + ".7z " + final_fold + " -m0=PPMD";
                    System.out.println(command);

                    String line;
                    Process p = Runtime.getRuntime().exec(command);

                    BufferedReader bri = new BufferedReader
                            (new InputStreamReader(p.getInputStream()));

                    BufferedReader bre = new BufferedReader
                            (new InputStreamReader(p.getErrorStream()));

                    while ((line = bri.readLine()) != null) {
                        System.out.println(line);
                    }
                    bri.close();

                    while ((line = bre.readLine()) != null) {
                        System.out.println(line);
                    }
                    bre.close();

                    p.waitFor();

                    Date endDate = new Date();
                    long taskCpuTimeNano = getCpuTime() - startCpuTimeNano;

                    double diff = endDate.getTime() - startDate.getTime();
                    double diffSeconds = diff / 1000;
                    double diffMinutes = diff / (60 * 1000);

                    System.out.println("Time in seconds*: " + (double) taskCpuTimeNano / 1000000000.0 + " seconds.");
                    System.out.println("Time in seconds: " + diffSeconds + " seconds.");
                    System.out.println("Time in minutes: " + diffMinutes + " minutes.");
                }
            }
        } 

    }
 
}
