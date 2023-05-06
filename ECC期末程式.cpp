#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int length=0;
char ch;
int main()
{
	FILE *receivedata = fopen("receive_data.txt", "r"); //讀取receive_data.txt檔
	if (receivedata == NULL){  //若無法讀取顯示error to read
		printf("error to read \n");
	}
	else
	{
		int receive[3000];  //預設data長度範圍
		printf("receive_data:\n");
		while((ch=getc(receivedata))!=EOF){   //從FILE receivedata取值
			int readdata=atoi(&ch);  //atoi將字串裡的數字字元轉化為整形數。返回整形值。
			receive[length] = readdata;
			printf("%d",receive[length]);
			length+=1;
		}
	    /*for(int i=0;i<length;i+=3){
			printf("%d%d%d ",receive[i],receive[i+1],receive[i+2]);
		}*/
		printf("\n");
		//printf("%d",length);
		//printf("\n");
		//times=length/3;
		int memory[7]={0};   //共有七個memory並先初始化為0
		int output[128][3];  //輸出為128(2的七次方)*3的矩陣型態,code rate為1/3 output的值表示經過generation matrix運算後的值
		int mem_state[128][7]={0};
		int x=0;
		//decode  建立table 找出所有的state
		for (int i=0; i<64; i++){    //i++，由state(i=1),...,state(i=63)各自生成兩個state，因此共是0~127個state。 (64*2=128)
			for(int input=0; input<2; input++){  //input為2binary,故只有0和1
				for(int y=0;y<7;y++){
					memory[y]=mem_state[i][y];
				}      //根據Optimum convolution codes with rate R=1/3, v=7時,g0=225 g1=331 g2=367 ,dfree=16 r為7.27(dB)
				output[x][0]=input^memory[2]^memory[4]^memory[6];  //g0=225  10010101--->246  (二進轉八進計算)
				output[x][1]=input^memory[0]^memory[2]^memory[3]^memory[6];  //g1=331  11011001--->0236
				output[x][2]=input^memory[0]^memory[1]^memory[2]^memory[4]^memory[5]^memory[6];  //g2=367  11110111--->012456
				memory[6]=memory[5];  //shift  input從state的第一個記憶體開始，一個一個shift
				memory[5]=memory[4];
				memory[4]=memory[3];
				memory[3]=memory[2];
				memory[2]=memory[1];
				memory[1]=memory[0];
				memory[0]=input;
				for (int y=0; y<7; y++){
					mem_state[x][y]=memory[y];
				}
				x++;
			}
		}

		int weight[1000][128][2];
		int k;

		//紀錄State與output的weight
		for (int i=0; i<(length/3); i++){  //sum weight  由於R=1/3, 將count/3 每三個bits去做weight 計算
			if(i>6){
				for(int x=0; x<128; x++){		//x的範圍最大為128
				k=3*i;
					for (int y=0; y<3; y++){
						if(receive[k]!=output[x][y]^1){
							weight[i][x][0]=weight[i][x][0]+1;
						}
						k++;
					}
					weight[i][x][1]=3-weight[i][x][0];    // 因為code rate=1/3,兩條codeword branch最大的誤差為3，故用3去減
				}
			}
			else{      //考慮 i<7 的情況
				for(int x=0; x<pow(2,i+1); x++){     //x的範圍最大為2的i+1次方 即 0<=x<2的(0~6)+1次方
					k=3*i;
					for (int y=0; y<3; y++){
						if(receive[k]!=output[x][y]){    //如果 receive的值和 output的值不一樣  weight會增加一 (output的值表示經由generation matrix運算後的值)
							weight[i][x][0]=weight[i][x][0]+1;
						}
						k++;
					}
				}
			}
		}
		/*for (int i=0; i<10; i++){
			for (int j=0; j<64; j++){
				printf("%d",weight[i][j][0]);
			}
			printf("\n");
		}*/

		int route[1000][128];  //路徑
		for (int i=1; i<(length/3); i++){
			//由前兩種weight相比
			if(i<7){     //每個state只有唯一一條路徑路徑
				for(int x=0; x<pow(2,i); x++){
					for (int y=2*x; y<2+2*x; y++){
						weight[i][y][0]=weight[i-1][x][0]+weight[i][y][0];
						//紀錄所選擇的route
						route[i][y]=x;  //此route的矩陣是指在i時，這個state是從上一個i的哪種狀態過來的
					}
				}
			}
			else{    //比較同一state不同來源的weight 並紀錄來自哪個state
				for(int x=0; x<128; x=x+2){
					for(int y=0; y<2; y++){
						weight[i][x+y][1]=weight[i][x+y][1]+weight[i-1][x/2][0];
						weight[i][x+y][0]=weight[i][x+y][0]+weight[i-1][(x/2)+64][0];
						if(weight[i][x+y][0] > weight[i][x+y][1]){   //比較各個weight
							weight[i][x+y][0]=weight[i][x+y][1];
							route[i][x+y]=x/2;
						}
						else{
							route[i][x+y]=64+x/2;
						}
					}
				}
			}
		}

		//printf("%d",length);
		//找出最佳route解出receive code
		int back;
		int decode[1000];

		for (int i=(length/3)-1; i>0; i--){  //從後面往前連接起來，找出最佳路徑
			//printf("%d",length);
			back=route[i][back];   //最後一個時間i-1時，是從state0開始找回去，其值為i-2時的state
			//printf("%d",length);
			if(mem_state[back][0]!=0){  //若 mem_state為0，解碼判斷為1
				decode[i-1]=1;
			}
			else{
				decode[i-1]=0;
			}
		}

		/*for(int i=0; i<(length/3); i++){
			printf("%d",decode[i]);
		}
		printf("\n");*/
		//encode計算

		int encode[3000][3];
		int memory2[7]={0};
		for(int i=0; i<(length/3); i++){     //注意此時的encode為state2非state 但仍使用相同的 g0,g1,g2
			encode[i][0]=decode[i]^memory2[2]^memory2[4]^memory2[6];
			encode[i][1]=decode[i]^memory2[0]^memory2[2]^memory2[3]^memory2[6];
			encode[i][2]=decode[i]^memory2[0]^memory2[1]^memory2[2]^memory2[4]^memory2[5]^memory2[6];
			memory2[6]=memory2[5];   //shift
			memory2[5]=memory2[4];
			memory2[4]=memory2[3];
			memory2[3]=memory2[2];
			memory2[2]=memory2[1];
			memory2[1]=memory2[0];
			memory2[0]=decode[i];
		}
		printf("encode_data:\n");
		for(int i=0; i<(length/3); i++){
			printf("%d%d%d ",encode[i][0],encode[i][1],encode[i][2]);
		}
		//計算error總數
		int error=0;     //用error數量記錯誤數量
		printf("\n");
		for(int i=0; i<(length/3); i++){
			for(int j=0; j<3; j++){
				if(receive[3*i+j]!=encode[i][j]){
					error = error+1;
					printf("error position: %d \n", 3*i+j+1);
				}
			}
		}
		printf("total error number= %d\n",error);
	}

	system("pause");
	return 0;
}
