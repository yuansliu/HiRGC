# include <cstdlib>
# include <iostream>
# include <fstream>
# include <ctime>
# include <sys/time.h>
# include <cmath>
# include <vector>
# include <algorithm>
# include <map>
# include <climits>
# include <cstdint>
# include <cstring>
# include <string>
# include <utility>
# include <unistd.h>

using namespace std;

const int MAX_CHAR_NUM = 1<<28;
// const unsigned char code_rule[4] = {0, 1, 2, 3};//A-0; C-1; G-2; T-3;
const char invert_code_rule[4] = {'A', 'C', 'G', 'T'};

struct POSITION_RANGE { //record lower case range
	int begin, length;//
};

struct POSITION_OTHER_CHAR {
	int pos, ch;
};

struct REF_NODE {
	int v, next;
};

//for target and ref genomic
char *meta_data;
POSITION_RANGE *pos_vec, *n_vec;
POSITION_OTHER_CHAR *other_char_vec;
int *line_break_vec;
char *ref_seq, *tar_seq;
int other_char_len, str_len, tar_seq_len, ref_seq_len, pos_vec_len, n_vec_len, line_break_len;

const int sub_str_num = 20;

void initial() {
	meta_data = new char[50];
	tar_seq = new char[MAX_CHAR_NUM];
	ref_seq = new char[MAX_CHAR_NUM];
}

void clear() {
	delete[] meta_data;
	delete[] tar_seq;
	delete[] ref_seq;
}

void readRefFile(char *refFile) {
	ref_seq_len = 0;
	char ch[1000];
	FILE *fp = fopen(refFile, "r");
	if (NULL == fp) {
		printf("fail to open file %s\n", refFile);
		return;
	}
	int temp_len;
	char temp_ch;
	while (fscanf(fp, "%s", ch) != EOF) {
		if (ch[0] == '>' || ch[0] == '@') {
			continue;
		}
		// printf("%s\n", ch);
		temp_len = strlen(ch);
		for (int i = 0; i < temp_len; i++) {
			temp_ch = ch[i];
			if (islower(temp_ch)) {
				temp_ch = toupper(temp_ch);
			}
			if (temp_ch == 'A' || temp_ch == 'G' || temp_ch == 'C' || temp_ch == 'T' ) {
				ref_seq[ref_seq_len++] = temp_ch;
			}
		}
	}
	ref_seq[ref_seq_len] = '\0';
	fclose(fp);
}

void readRunLengthCoding(FILE *fp, int &vec_len, int **vec) {
	int code_len, temp;
	vector<int> code;
	fscanf(fp, "%d", &code_len);
	for (int i = 0; i < code_len; i++) {
		fscanf(fp, "%d", &temp);
		code.push_back(temp);
	}

	vec_len = 0;
	for (int i = 1; i < code_len; i += 2) {
		vec_len += code[i];
	}

	if (vec_len > 0) {
		*vec = new int[vec_len];
		int k = 0;
		for (int i = 0; i < code_len; i += 2) {
			for (int j = 0; j < code[i+1]; j++) {
				(*vec)[k++] = code[i];
			}
		}
	}
}

void readOtherData(FILE *fp) {
	meta_data = new char[100];
	fscanf(fp, "%s", meta_data);

	readRunLengthCoding(fp, line_break_len, &line_break_vec);

	fscanf(fp, "%d", &pos_vec_len);
	// printf("pos_vec_len: %d\n", pos_vec_len);
	if (pos_vec_len > 0) {
		pos_vec = new POSITION_RANGE[pos_vec_len];
		for (int i = 0; i < pos_vec_len; i++) {
			fscanf(fp, "%d%d", &pos_vec[i].begin, &pos_vec[i].length);
		}
	}

	fscanf(fp, "%d", &n_vec_len);
	// printf("n_vec_len: %d\n", n_vec_len);
	if (n_vec_len > 0) {
		n_vec = new POSITION_RANGE[n_vec_len];
		for (int i = 0; i < n_vec_len; i++) {
			fscanf(fp, "%d%d", &n_vec[i].begin, &n_vec[i].length);
		}
	}

	fscanf(fp, "%d", &other_char_len);
	// printf("other_char_len: %d\n", other_char_len);
	if (other_char_len > 0) {
		other_char_vec = new POSITION_OTHER_CHAR[other_char_len];
		for (int i = 0; i < other_char_len; i++) {
			fscanf(fp, "%d", &other_char_vec[i].pos);
		}
		//read other_char_vec information
		vector<int> arr;
		int size, temp;
		fscanf(fp, "%d", &size);
		for (int i = 0; i < size; i++) {
			fscanf(fp, "%d", &temp);
			arr.push_back(temp);
		}
		if (size < 10) {
			char *temp_str = new char[other_char_len+10];
			fscanf(fp, "%s", temp_str);
			for (int i = 0; i < other_char_len; i++) {
				other_char_vec[i].ch = arr[temp_str[i]-'0'];
			}
			delete[] temp_str;
		} else {
			int bit_num = ceil(log(size)/log(2));
			int v_num = floor(32.0/bit_num);

			for (int i = 0; i < other_char_len; ) {
				unsigned int v;
				fscanf(fp, "%u", &v);
				vector<int> temp_arr;

				int temp_i = i;
				for (int j = 0; j < v_num && temp_i < other_char_len; j++, temp_i++) {
					int mod = v%(1<<bit_num);
					v >>= bit_num;
					temp_arr.push_back(arr[mod]);
				}
				for (int j = temp_arr.size()-1; j >= 0 && i < other_char_len; j--, i++) {
					other_char_vec[i].ch = temp_arr[j];
				}
			}
		}
	}
}

bool exitSpace(char *str) {
	for (int i = 0; str[i]; i++) {
		if (str[i] == ' ') {
			return true;
		}
	}
	return false;
} 

void readCompressedFile(char *tarFile) {
	FILE *fp = fopen(tarFile, "r");
	if (NULL == fp) {
		printf("ERROR! fail to open file %s\n", tarFile);
	}
	readOtherData(fp);

	int pre_pos = 0, cur_pos, length;;
	tar_seq_len = 0;
	// printf("%s\n", ref_seq);
	char *str = new char[20000];
	while (fgets(str, MAX_CHAR_NUM, fp) != NULL) {
		if (str[0] == '\n') continue;
		// printf("cmd: %s", str);
		if (exitSpace(str)) {
			sscanf(str, "%d%d", &cur_pos, &length);
			cur_pos += pre_pos;
			length += sub_str_num;
			pre_pos = cur_pos + length;
			for (int i = cur_pos, j = 0; j < length; j++, i++) {
				tar_seq[tar_seq_len++] = ref_seq[i];
			}
			tar_seq[tar_seq_len] = '\0';
		} else {
			int str_len = strlen(str);
			for (int i = 0; i < str_len-1; i++) {
				tar_seq[tar_seq_len++] = invert_code_rule[str[i]-'0'];
			}
		}
	}
	tar_seq[tar_seq_len] = '\0';
	// printf("%s\n", tar_seq);
	fclose(fp);
	delete[] str;
}
 
void saveDecompressedData(char *resultFile) {
	FILE *fp = fopen(resultFile, "w");
	if (NULL == fp) {
		printf("ERROR! fail to open file %s\n", resultFile);
		return;
	}

	for (int i = 1; i < other_char_len; i++) {
		other_char_vec[i].pos += other_char_vec[i-1].pos;
	}

	char *temp_seq = new char[MAX_CHAR_NUM];
	
	strcpy(temp_seq, tar_seq);
	int tt = 0, j = 0;

	for (int i = 0; i < other_char_len; i++) {
		while (tt < other_char_vec[i].pos && tt < tar_seq_len) {
			tar_seq[j++] = temp_seq[tt++];
		}
		tar_seq[j++] = other_char_vec[i].ch + 'A';
	}
	while (tt < tar_seq_len) {
		tar_seq[j++] = temp_seq[tt++];
	}
	tar_seq[j] = '\0';
	tar_seq_len = j;

	int str_len = 0;
	int r = 0; 

	char *str = new char[MAX_CHAR_NUM];

	for (int i = 0; i < n_vec_len; i++) {
		for (int j = 0; j < n_vec[i].begin; j++) {
			str[str_len++] = tar_seq[r++];
		}
		for (int j = 0; j < n_vec[i].length; j++) {
			str[str_len++] = 'N';
		}
	}
	while (r < tar_seq_len) {
		str[str_len++] = tar_seq[r++];
	}
	str[str_len] = '\0';

	fprintf(fp, "%s\n", meta_data);

	int k = 0;
	for (int i = 0; i < pos_vec_len; i++) {
		k += pos_vec[i].begin;
		int temp = pos_vec[i].length;
		for (int j = 0; j < temp; j++) {
			str[k] = tolower(str[k]);
			k++;
		}
	}
	
	for (int i = 1; i < line_break_len; i++) {
		line_break_vec[i] += line_break_vec[i-1];
	}
	int k_lb = 0;
	
	int temp_seq_len = 0;
	for (int i = 0; i < str_len; i++) {
		if(i == line_break_vec[k_lb]) {
			temp_seq[temp_seq_len++] = '\n';
			k_lb++;
		}
		temp_seq[temp_seq_len++] = str[i];
	}
	temp_seq[temp_seq_len] = '\0';
	fprintf(fp, "%s", temp_seq);
	// printf("k_lb: %d\n", k_lb);
	while ((k_lb++) < line_break_len) {
		fprintf(fp, "\n");
	}
	fclose(fp);

	delete[] temp_seq;
	delete[] str;

	if (pos_vec_len > 0) delete[] pos_vec;
	delete[] line_break_vec;
	if (n_vec_len > 0) delete[] n_vec;
	if (other_char_len > 0) delete[] other_char_vec;
}

void decompressFile(char *refFile, char *tarFile, char *resultFile) {
	readRefFile(refFile);//must read first;
	readCompressedFile(tarFile);
	saveDecompressedData(resultFile);
}

void decompressFile(char *refFile, char *tarFile) {
	initial();//****important
	char cmd[200]; 
	sprintf(cmd, "./7za x %s", tarFile);
	system(cmd);

	printf("decompressing...\n");
	
	char temp[50];
	sprintf(temp, "%s", tarFile);
	for (int i = 0; temp[i]; i++) {
		if (temp[i] == '.' && temp[i+1] == '7') {
			temp[i] = '\0';
			break;
		}
	}
	char res[100];
	sprintf(res, "dec_%s", temp);
	decompressFile(refFile, temp, res);

	clear();
}

void decompressGenome(char *ref_fold, char *tar_7z, vector<string> &chr_name_list) {
	initial();//****important
	char ref[100], tar[100], res[100];
	char cmd[200]; 
	sprintf(cmd, "./7za x %s", tar_7z);
	system(cmd);

	char temp[10];
	sprintf(temp, "%s", tar_7z);
	for (int i = 0; temp[i]; i++) {
		if (temp[i] == '.' && temp[i+1] == '7') {
			temp[i] = '\0';
			break;
		}
	}

	// sprintf(cmd, "rm -rf %s_dec", temp);
	// system(cmd);

	sprintf(cmd, "mkdir %s_dec", temp);
	system(cmd);
	printf("the decompress result saved in '%s_dec'\n", temp);
	
	int size = chr_name_list.size();
	for (int i = 0; i < size; i++) {
		sprintf(ref, "%s/%s", ref_fold, chr_name_list[i].c_str());
		sprintf(tar, "%s/%s", temp, chr_name_list[i].c_str());
		printf("decompressing %s\n", tar);
		sprintf(res, "%s_dec/%s", temp, chr_name_list[i].c_str());

		decompressFile(ref, tar, res);
	}
	
	sprintf(cmd, "rm -rf %s", temp);
	system(cmd);

	clear();
}

vector<string> defalt_name_list = {"chr1.fa", "chr2.fa", "chr3.fa", "chr4.fa", 
                    "chr5.fa", "chr6.fa", "chr7.fa", "chr8.fa", "chr9.fa", "chr10.fa", 
                    "chr11.fa", "chr12.fa", "chr13.fa", "chr14.fa", "chr15.fa", "chr16.fa", "chr17.fa", 
                    "chr18.fa", "chr19.fa", "chr20.fa", "chr21.fa", "chr22.fa", "chrX.fa", "chrY.fa"};

/*int main() {
	char refFile[100], tarFile[100], resultFile[100];
	initial();
	sprintf(refFile, "HG18/chr1.fa");
	sprintf(tarFile, "temp.cd");
	sprintf(resultFile, "de.fa");
	decompressFile(refFile, tarFile, resultFile);
	return 0;
}*/

/*int main(int argc, char* argv[]) {
	//1 mapping_file, 2 ref_fold, 3 7z_file
	char ref[100], tar[100], res[100], str1[100], str2[100];
	char cmd[200]; 
	if (argc != 4) {
		printf("parameters not correct!\n");
	} else {
		FILE *fp_map = fopen(argv[1], "r");
		if (!fp_map) {
			printf("fail to open mapping file.");
		} else {
			initial();//important

			sprintf(cmd, "./7za x %s", argv[3]);
			system(cmd);

			char temp[10];
			sprintf(temp, "%s", argv[3]);
			for (int i = 0; temp[i]; i++) {
				if (temp[i] == '.' && temp[i+1] == '7') {
					temp[i] = '\0';
					break;
				}
			}

			sprintf(cmd, "rm %s_dec", temp);
			system(cmd);

			sprintf(cmd, "mkdir %s_dec", temp);
			system(cmd);
			printf("the decompress result saved in '%s_dec'\n", temp);
			
			struct  timeval  start;
			struct  timeval  end;
			unsigned long timer;
			gettimeofday(&start,NULL);
			
			while(fscanf(fp_map, "%s%s", str1, str2) != EOF) {
				sprintf(ref, "%s/%s", argv[2], str1);
				sprintf(tar, "%s/%s", temp, str2);
				printf("decompressing %s\n", tar);
				sprintf(res, "%s_dec/%s", temp, str2);

				decompressFile(ref, tar, res);
			}
			
			sprintf(cmd, "rm -rf %s", temp);
			system(cmd);

			gettimeofday(&end,NULL);
			timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
			printf("total compression timer = %lf ms\n", timer/1000.0);
		}
	}
}*/

bool getList(char *list_file, vector<string> &name_list) {
	FILE *fp = fopen(list_file, "r");
	if (fp == NULL) {
		printf("%s open fail!\n", list_file);
		return false;
	}
	char str[100];
	while (fscanf(fp, "%s", str) != EOF) {
		name_list.push_back(string(str));
	}
	fclose(fp);
	if (name_list.size() == 0) {
		printf("%s is empty!\n", list_file);
		return false;
	}
	return true;
}

void setDefaltName(vector<string> &name_list) {
	int size = defalt_name_list.size();
	for (int i = 0; i < size; i++) {
		name_list.push_back(defalt_name_list[i]);
	}
}

int main(int argc, char *argv[]) {
	vector<string> chr_name_list;

	bool flag = true, decompressed = false;
	int oc;
	char *mode = NULL, *ref_file = NULL, *tar_file = NULL, *ref_fold = NULL, *tar_fold = NULL;

	struct  timeval  start;
	struct  timeval  end;
	unsigned long timer;
	gettimeofday(&start,NULL);

	if ((oc = getopt(argc, argv, "m:")) >= 0) {
		mode = optarg;
	}

	if (strcmp(mode, "file") == 0) {
		while ((oc = getopt(argc, argv, "r:t:")) >= 0) {
			switch(oc) {
				case 'r':
					ref_file = optarg;
					printf("%s\n", ref_file);
					break;
				case 't':
					tar_file = optarg;
					break;
				case '?':
					flag = false;
					break;
			}
		}

		if (flag && ref_file && tar_file) {
			decompressFile(ref_file, tar_file);
			decompressed = true;
		}

	} else 
	if (strcmp(mode, "genome") == 0 ){
		while ((oc = getopt(argc, argv, "r:t:n:")) >= 0) {
			switch(oc) {
				case 'r':
					ref_fold = optarg;
					break;
				case 't':
					tar_fold = optarg;
					break;
				case 'n':
					if (strcmp(optarg, "default")) {
						setDefaltName(chr_name_list);
					} else {
						flag &= getList(optarg, chr_name_list);
					}
					break;
				case '?':
					flag = false;
					break;
			}
		}

		if (flag && ref_fold && tar_fold && chr_name_list.size() > 0) {
			decompressGenome(ref_fold, tar_fold, chr_name_list);
			decompressed = true;
		}

	} 

	if (!decompressed) {
		printf("arguments error...\n");
		return 0;
	}

	gettimeofday(&end,NULL);
	timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
	printf("total decompression timer = %lf ms; %lf min\n", timer/1000.0, timer/1000.0/1000.0/60.0);
}

