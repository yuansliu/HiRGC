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
const int code_rule[4] = {0, 1, 2, 3};//A-0; C-1; G-2; T-3;
const char invert_code_rule[4] = {'A', 'C', 'G', 'T'};
const int max_arr_num_bit = 30;
const int max_arr_num_bit_shift = max_arr_num_bit>>1;
const int max_arr_num = 1<<max_arr_num_bit;
const int min_size = 1<<20;

const int sub_str_num = 20;

struct POSITION_RANGE { //record lower case range
	int begin, length;//
};

struct POSITION_OTHER_CHAR {
	int pos, ch;
};

//for target and ref genomic
char *meta_data;
POSITION_RANGE *pos_vec, *n_vec;
POSITION_OTHER_CHAR *other_char_vec;
int *line_break_vec;
int *tar_seq_code, *ref_seq_code;
int other_char_len, str_len, tar_seq_len, ref_seq_len, pos_vec_len, n_vec_len, line_break_len;
int *loc;
int *point;
char *dismatched_str;

inline void initial() {
	meta_data = new char[50];
	pos_vec = new POSITION_RANGE[min_size];
	line_break_vec = new int[1<<23];
	n_vec = new  POSITION_RANGE[min_size];
	other_char_vec = new POSITION_OTHER_CHAR[min_size];

	tar_seq_code = new int[MAX_CHAR_NUM];
	ref_seq_code = new int[MAX_CHAR_NUM];

	loc = new int[MAX_CHAR_NUM];
	point = new int[max_arr_num];

	dismatched_str = new char[min_size];
}

inline void clear() {
	delete[] meta_data;
	delete[] pos_vec;
	delete[] line_break_vec;
	delete[] n_vec;
	delete[] other_char_vec;
	delete[] tar_seq_code;
	delete[] ref_seq_code;
	delete[] loc;
	delete[] point;
	delete[] dismatched_str;
}

int agctIndex(char ch) {
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
    return 4;
}

void readRefFile(char *refFile) {
	int _ref_seq_len = 0;
	char ch[128];
	FILE *fp = fopen(refFile, "r");
	if (NULL == fp) {
		printf("fail to open file %s\n", refFile);
		return;
	}
	int temp_len, index;
	char temp_ch;
	
	fscanf(fp, "%s", ch);//meta_data

	while (fscanf(fp, "%s", ch) != EOF) {
		temp_len = strlen(ch);
		for (int i = 0; i < temp_len; i++) {
			temp_ch = ch[i];
			if (islower(temp_ch)) {
				temp_ch = toupper(temp_ch);
			}
			index = agctIndex(temp_ch);
			if (index^4) {
				ref_seq_code[_ref_seq_len++] = index;
			}
		}
	}
	fclose(fp);
	ref_seq_len = _ref_seq_len;
}

void readTarFile(char *tarFile) {
	FILE *fp = fopen(tarFile, "r");
	if (NULL == fp) {
		printf("fail to open file %s\n", tarFile);
		return;
	}
	char ch[128], chr;
	
	int _tar_seq_len = 0;
	other_char_len = 0, n_vec_len = 0;
	pos_vec_len = line_break_len = 0;

	int letters_len = 0, n_letters_len = 0, index, ch_len;
	bool flag = true, n_flag = false; // first is upper case //first is not N

	fscanf(fp, "%s", meta_data); //meta_data

	while (fscanf(fp, "%s", ch) != EOF) {
		ch_len = strlen(ch);

		for (int i = 0; i < ch_len; i++) {
			chr = ch[i];
			if (islower(chr)) {
				if (flag) { //previous is upper case
					flag = false; //change status of switch
					pos_vec[pos_vec_len].begin = (letters_len);
					letters_len = 0;
				}
				chr = toupper(chr);
			} else {
				if (isupper(chr)) {
					if (!flag) {
						flag = true;
						pos_vec[pos_vec_len].length = (letters_len);
						pos_vec_len++;
						letters_len = 0;
					}
				}
			}
			letters_len++;

			//ch is an upper letter
			if (chr != 'N') {
				index = agctIndex(chr);
				if (index^4) {
					// tar_seq_code[tar_seq_len++] = code_rule[index];
					tar_seq_code[_tar_seq_len++] = index;
				} else {
					other_char_vec[other_char_len].pos = _tar_seq_len;
					other_char_vec[other_char_len].ch = chr-'A';
					other_char_len++;
				}
			}

			if (!n_flag) {
				if (chr == 'N') {
					n_vec[n_vec_len].begin = n_letters_len;
					n_letters_len = 0;
					n_flag = true;
				}
			} else {//n_flag = true
				if (chr != 'N'){
					n_vec[n_vec_len].length = n_letters_len;
					n_vec_len++;
					n_letters_len = 0;
					n_flag = false;
				}
			}
			n_letters_len++;

		}
		line_break_vec[line_break_len++] = ch_len;
	}

	if (!flag) {
		pos_vec[pos_vec_len].length = (letters_len);
		pos_vec_len++;
	}
	
	if (n_flag) {
		n_vec[n_vec_len].length = (n_letters_len);
		n_vec_len++;
	}

	for (int i = other_char_len-1; i > 0; i--) {
		other_char_vec[i].pos -= other_char_vec[i-1].pos;
	}
	fclose(fp);
	tar_seq_len = _tar_seq_len;
}

void writeRunLengthCoding(FILE *fp, int vec_len, int *vec) {
	vector<int> code;
	if (vec_len > 0) {
		code.push_back(vec[0]);
		int pre_value = code[0], cnt = 1;
		for (int i = 1; i < vec_len; i++) {
			if (vec[i] == pre_value) {
				cnt++;
			} else {
				code.push_back(cnt);
				code.push_back(vec[i]);
				pre_value = vec[i];
				cnt = 1;
			}
		}
		code.push_back(cnt);
	}
	int code_len = code.size();
	fprintf(fp, "%d", code_len);

	for (int i = 0; i < code_len; i++) {
		fprintf(fp, " %d", code[i]);
	}
	fprintf(fp, "\n");
}

void saveOtherData(FILE *fp) {
	fprintf(fp, "%s\n", meta_data);
	//------------------------
	writeRunLengthCoding(fp, line_break_len, line_break_vec);
	//------------------------
	fprintf(fp, "%d", pos_vec_len);
	for (int i = 0; i < pos_vec_len; i++) {
		fprintf(fp, " %d %d", pos_vec[i].begin, pos_vec[i].length);
	}
	//------------------------
	fprintf(fp, "\n%d", n_vec_len);
	for (int i = 0; i < n_vec_len; i++) {
		fprintf(fp, " %d %d", n_vec[i].begin, n_vec[i].length);
	}
	//------------------------
	fprintf(fp, "\n%d", other_char_len);
	if (other_char_len > 0) {
		int flag[30];
		for (int i = 0; i < 26; i++) {
			flag[i] = -1;
		}
		vector<int> arr;
	
		for (int i = 0; i < other_char_len; i++) {
			fprintf(fp, " %d", other_char_vec[i].pos);

			int temp = other_char_vec[i].ch;
			if (flag[temp] == -1) {
				arr.push_back(temp);
				flag[temp] = arr.size()-1;
			}
		}

		//save other char information
		int size = arr.size();
		fprintf(fp, " %d", size);
		for (int i = 0; i < size; i++) {
			fprintf(fp, " %d", arr[i]);
		}

		if (size < 10) {
			fprintf(fp, " ");
			for (int i = 0; i < other_char_len; i++) {
				fprintf(fp, "%d", flag[other_char_vec[i].ch]);
			}
		} else {
			int bit_num = ceil(log(size)/log(2));
			int v_num = floor(32.0/bit_num);

			for (int i = 0; i < other_char_len; ) {
				// fprintf(fp, "%d ", flag[other_char_vec[i].ch]);
				unsigned int v = 0;
				for (int j = 0; j < v_num && i < other_char_len; j++, i++) {
					v <<= bit_num;
					v += flag[other_char_vec[i].ch];
				}
				fprintf(fp, " %u", v);
			}
		}
	}
	fprintf(fp, "\n");
}

void preProcessRef() {
	for (int i = 0; i < max_arr_num; i++) {
		point[i] = -1;
	}

	// memset(point, -1, sizeof(int)*max_arr_num);

	uint64_t value = 0;
	for (int k = sub_str_num - 1; k >= 0; k--) {
		value <<= 2;
		value += ref_seq_code[k];
	}
	int id = value&(uint64_t)(max_arr_num-1);
	loc[0] = point[id];
	point[id] = 0;

	int step_len = ref_seq_len - sub_str_num + 1;
	int shift_bit_num = (sub_str_num*2-2);
	int one_sub_str = sub_str_num - 1;

	for (int i = 1; i < step_len; i++) {
		value >>= 2;
		value += ((uint64_t)ref_seq_code[i + one_sub_str]<<shift_bit_num);
		
		id = value&(uint64_t)(max_arr_num-1);
		loc[i] = point[id];
		point[id] = i;
	}	
}

void searchMatch(char *refFile) {
	FILE *fp = fopen(refFile, "w");
	if (NULL == fp) {
		printf("ERROR! fail to open file %s\n", refFile);
		return;
	}

	saveOtherData(fp);

	int pre_dismatch_point = 0, pre_pos = 0;
	int step_len = tar_seq_len - sub_str_num + 1;
	int max_length, max_k;

	int dis_str_len = 0, i, id, k, ref_idx, tar_idx, length, cur_pos;;
	uint64_t tar_value;

	for (i = 0; i < step_len; i++) {
		tar_value = 0;
		for (k = sub_str_num - 1; k >= 0; k--) {
			tar_value <<= 2;
			tar_value += tar_seq_code[i+k];
		}

		id = point[tar_value&(uint64_t)(max_arr_num-1)];
		if (id > -1) {
			max_length = -1;
			max_k = -1;
			for (k = id; k != -1; k = loc[k]) {
					//[pos[k], pos[k]+sub_str_num-1] //at least max_arr_num_bit_shift
				ref_idx = k + max_arr_num_bit_shift;//loc[k].pos == k
				tar_idx = i + max_arr_num_bit_shift;
				length = max_arr_num_bit_shift;
				while (ref_idx < ref_seq_len && tar_idx < tar_seq_len && ref_seq_code[ref_idx++] == tar_seq_code[tar_idx++]) {
					length++;
				}
				if (length >= sub_str_num && length > max_length) {
					max_length = length;
					max_k = k;
				}
			}
			if (max_k > -1) {
				//first print mismatch substring
				if (dis_str_len > 0) {
					dismatched_str[dis_str_len] = '\0';
					fprintf(fp, "%s\n", dismatched_str);
					dis_str_len = 0;
				}
				//then print match substring
				cur_pos = max_k - pre_pos;
				pre_pos = max_k + max_length;

				fprintf(fp, "%d %d\n", cur_pos, max_length - sub_str_num);

				i += max_length;
				pre_dismatch_point = i;
				if (i < tar_seq_len) {
					dismatched_str[dis_str_len++] = '0'+tar_seq_code[i];
				}
				continue;
			}
		}
		dismatched_str[dis_str_len++] = '0'+tar_seq_code[i];
	}

	for(; i < tar_seq_len; i++) {
		dismatched_str[dis_str_len++] = '0'+tar_seq_code[i];
	}
	if (dis_str_len > 0) {
		dismatched_str[dis_str_len] = '\0';
		fprintf(fp, "%s\n", dismatched_str);
	}
	fclose(fp);
	// gettimeofday(&end,NULL);
	// timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
	// printf("compress timer = %lf ms\n", timer/1000.0);
}

inline void compressFile(char *refFile, char *tarFile, char *resFile) {
	readRefFile(refFile);
	preProcessRef();
	readTarFile(tarFile);
	searchMatch(resFile);
}

/*int main() {
	sub_str_num = 19;
	char refFile[100], tarFile[100], resFile[100];
	initial();
	sprintf(refFile, "HG18/chr1.fa");
	sprintf(tarFile, "HG19/chr1.fa");	
	// sprintf(refFile, "simu-ko-131.fa");
	// sprintf(tarFile, "simu-ko-224.fa");
	sprintf(resFile, "temp.cd");
	compressFile(refFile, tarFile, resFile);
	return 0;
}*/

/*int main(int argc, char* argv[]) {//hirgc mapping.txt refFold tarFold
	char ref[100], tar[100], res[100], str1[100], str2[100];
	char cmd[200]; 

	if (argc != 4) {
		printf("parameters not correct!\n");
	} else {
		FILE *fp_map = fopen(argv[1], "r");
		if (!fp_map) {
			printf("fail to open mapping file.");
		} else {
			sub_str_num = 20; //important*************
			initial(); //important*************

			char temp[10];
			sprintf(temp, "%s_ref_%s", argv[3], argv[2]);
			sprintf(cmd, "mkdir %s", temp);
			system(cmd);
			
			printf("sub_str_num: %d\n", sub_str_num);
			printf("max_arr_num_bit: %d\n", max_arr_num_bit);

			struct  timeval  start;
			struct  timeval  end;
			unsigned long timer;
			gettimeofday(&start,NULL);

			while(fscanf(fp_map, "%s%s", str1, str2) != EOF) {
				sprintf(ref, "%s/%s", argv[2], str1);
				sprintf(tar, "%s/%s", argv[3], str2);
				printf("compressing %s\n", tar);
				sprintf(res, "./%s/%s", temp, str2);

				compressFile(ref, tar, res);
			}
			
			sprintf(cmd, "./7za a %s.7z %s -m0=PPMd", temp, temp);
			system(cmd);
			// system("rm -r temp");
			gettimeofday(&end,NULL);
			timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
			printf("total compression timer = %lf ms\n", timer/1000.0);
		}
	}
	clear();
	return 0;
}*/

void compressFile(char *refFile, char *tarFile) {
	char res[100];
	char cmd[200];
	initial(); //important************* must initial before other operation

	sprintf(res, "%s_ref_%s", tarFile, refFile);
	compressFile(refFile, tarFile, res);
	sprintf(cmd, "./7za a %s.7z %s -m0=PPMd", res, res);
	system(cmd);
	sprintf(cmd, "rm %s", res);
	system(cmd);

	clear();
}

void compressGenome(char *refFold, char *tarFold, vector<string> &chr_name_list) {
    char ref[100], tar[100], res[100];
	char cmd[200];
	char temp[100];
	sprintf(temp, "%s_ref_%s", tarFold, refFold);
	sprintf(cmd, "rm -r %s", temp);
	system(cmd);
	sprintf(cmd, "mkdir %s", temp);
	system(cmd);

	initial(); //important************* must initial before other operation

	int size = chr_name_list.size();
	for (int i = 0; i < size; i++) {
		sprintf(ref, "%s/%s", refFold, chr_name_list[i].c_str());
		sprintf(tar, "%s/%s", tarFold, chr_name_list[i].c_str());
		printf("compressing %s\n", tar);
		sprintf(res, "%s/%s", temp, chr_name_list[i].c_str());
		compressFile(ref, tar, res);
	}
	
	sprintf(cmd, "./7za a %s.7z %s -m0=PPMd", temp, temp);
	system(cmd);

	sprintf(cmd, "rm -r %s", temp);
	system(cmd);
	
	clear();
}

int compressSet(char *ref_fold, vector<string> &fold_list, vector<string> &chr_name_list) {//hirgc refFold tarFold
	char ref[100], tar[100], res[100];
	char cmd[200]; 
	char temp[100];
	
	int fold_size = fold_list.size();

	for (int ti = 0; ti < fold_size; ti++) {
		sprintf(temp, "%s_ref_%s", fold_list[ti].c_str(), ref_fold);
		sprintf(cmd, "rm -r %s", temp);
		system(cmd);
		sprintf(cmd, "mkdir %s", temp);
		system(cmd);
	}

	initial(); //important*************

	int size = chr_name_list.size();
	for (int i = 0; i < size; i++) {

		sprintf(ref, "%s/%s", ref_fold, chr_name_list[i].c_str());

		readRefFile(ref);
		preProcessRef();

		for (int ti = 0; ti < fold_size; ti++) {
			sprintf(temp, "%s_ref_%s", fold_list[ti].c_str(), ref_fold);

			sprintf(tar, "%s/%s", fold_list[ti].c_str(), chr_name_list[i].c_str());
			readTarFile(tar);

			printf("compressing %s ...\n", tar);

			sprintf(res, "%s/%s", temp, chr_name_list[i].c_str());
			// compressFile(ref, tar, res);
			searchMatch(res);
		}
	}

	for (int ti = 0; ti < fold_size; ti++) {
		sprintf(temp, "%s_ref_%s", fold_list[ti].c_str(), ref_fold);

		sprintf(cmd, "./7za a %s.7z %s -m0=PPMd", temp, temp);
		system(cmd);

		sprintf(cmd, "rm -r %s", temp);
		system(cmd);
	}
	clear();
	return 0;
}

vector<string> defalt_name_list = {"chr1.fa", "chr2.fa", "chr3.fa", "chr4.fa", 
                    "chr5.fa", "chr6.fa", "chr7.fa", "chr8.fa", "chr9.fa", "chr10.fa", 
                    "chr11.fa", "chr12.fa", "chr13.fa", "chr14.fa", "chr15.fa", "chr16.fa", "chr17.fa", 
                    "chr18.fa", "chr19.fa", "chr20.fa", "chr21.fa", "chr22.fa", "chrX.fa", "chrY.fa"};

/*int main(int argc, char* argv[]) {//hirgc refFold tarFold
	char ref[100], tar[100], res[100];
	char cmd[200]; 

	string fold_list[] = {"YH", "HG18", "HG19", "HG38", "KOREF20090131", "KOREF20090224"};
	string name_list[] = {"chr1.fa", "chr2.fa", "chr3.fa", "chr4.fa", 
                    "chr5.fa", "chr6.fa", "chr7.fa", "chr8.fa", "chr9.fa", "chr10.fa", 
                    "chr11.fa", "chr12.fa", "chr13.fa", "chr14.fa", "chr15.fa", "chr16.fa", "chr17.fa", 
                    "chr18.fa", "chr19.fa", "chr20.fa", "chr21.fa", "chr22.fa", "chrX.fa", "chrY.fa"};
	// sub_str_num = 20; //important*************
	
	for (int ri = 0; ri < 1; ri++) {
		for (int ti = 0; ti < 6; ti++) {
			if (ri != ti) {
				char temp[100];
				sprintf(temp, "%s_ref_%s_test2", fold_list[ti].c_str(), fold_list[ri].c_str());
				sprintf(cmd, "mkdir %s", temp);
				system(cmd);

				struct  timeval  start;
				struct  timeval  end;
				unsigned long timer;
				gettimeofday(&start,NULL);

				initial(); //important*************

				for (int i = 0; i < 24; i++) {
					sprintf(ref, "%s/%s", fold_list[ri].c_str(), name_list[i].c_str());
					sprintf(tar, "%s/%s", fold_list[ti].c_str(), name_list[i].c_str());
					// printf("compressing %s\n", tar);
					sprintf(res, "%s/%s", temp, name_list[i].c_str());
					compressFile(ref, tar, res);
				}
				
				sprintf(cmd, "./7za a %s.7z %s -m0=PPMd", temp, temp);
				system(cmd);
				// system("rm -r temp");
				clear();
				
				gettimeofday(&end,NULL);
				timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
				// fprintf(out, "%s.7z: total compression timer = %lf ms; %lf min\n", temp, timer/1000.0, timer/1000.0/1000.0/60.0);
				printf("%s.7z: total compression timer = %lf ms; %lf min\n", temp, timer/1000.0, timer/1000.0/1000.0/60.0);

			}
		}
	}
	
	return 0;
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
	vector<string> fold_list;

	bool flag = true, compressed = false;
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
			compressFile(ref_file, tar_file);
			compressed = true;
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
			compressGenome(ref_fold, tar_fold, chr_name_list);
			compressed = true;
		}

	} else 
	if (strcmp(mode, "set") == 0) {
		while ((oc = getopt(argc, argv, "r:t:n:")) >= 0) {
			switch(oc) {
				case 'r':
					ref_fold = optarg;
					break;
				case 't':
					flag &= getList(optarg, fold_list);
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

		if (flag && ref_fold && fold_list.size() > 0 && chr_name_list.size() > 0) {
			compressSet(ref_fold, fold_list, chr_name_list);
			compressed = true;
		}
	}

	if (!compressed) {
		printf("arguments error...\n");
		return 0;
	}

	gettimeofday(&end,NULL);
	timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
	printf("total compression timer = %lf ms; %lf min\n", timer/1000.0, timer/1000.0/1000.0/60.0);
}
