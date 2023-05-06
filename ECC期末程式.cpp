#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int length=0;
char ch;
int main()
{
	FILE *receivedata = fopen("receive_data.txt", "r"); //Ū��receive_data.txt��
	if (receivedata == NULL){  //�Y�L�kŪ�����error to read
		printf("error to read \n");
	}
	else
	{
		int receive[3000];  //�w�]data���׽d��
		printf("receive_data:\n");
		while((ch=getc(receivedata))!=EOF){   //�qFILE receivedata����
			int readdata=atoi(&ch);  //atoi�N�r��̪��Ʀr�r����Ƭ���μơC��^��έȡC
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
		int memory[7]={0};   //�@���C��memory�å���l�Ƭ�0
		int output[128][3];  //��X��128(2���C����)*3���x�}���A,code rate��1/3 output���Ȫ�ܸg�Lgeneration matrix�B��᪺��
		int mem_state[128][7]={0};
		int x=0;
		//decode  �إ�table ��X�Ҧ���state
		for (int i=0; i<64; i++){    //i++�A��state(i=1),...,state(i=63)�U�ۥͦ����state�A�]���@�O0~127��state�C (64*2=128)
			for(int input=0; input<2; input++){  //input��2binary,�G�u��0�M1
				for(int y=0;y<7;y++){
					memory[y]=mem_state[i][y];
				}      //�ھ�Optimum convolution codes with rate R=1/3, v=7��,g0=225 g1=331 g2=367 ,dfree=16 r��7.27(dB)
				output[x][0]=input^memory[2]^memory[4]^memory[6];  //g0=225  10010101--->246  (�G�i��K�i�p��)
				output[x][1]=input^memory[0]^memory[2]^memory[3]^memory[6];  //g1=331  11011001--->0236
				output[x][2]=input^memory[0]^memory[1]^memory[2]^memory[4]^memory[5]^memory[6];  //g2=367  11110111--->012456
				memory[6]=memory[5];  //shift  input�qstate���Ĥ@�ӰO����}�l�A�@�Ӥ@��shift
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

		//����State�Poutput��weight
		for (int i=0; i<(length/3); i++){  //sum weight  �ѩ�R=1/3, �Ncount/3 �C�T��bits�h��weight �p��
			if(i>6){
				for(int x=0; x<128; x++){		//x���d��̤j��128
				k=3*i;
					for (int y=0; y<3; y++){
						if(receive[k]!=output[x][y]^1){
							weight[i][x][0]=weight[i][x][0]+1;
						}
						k++;
					}
					weight[i][x][1]=3-weight[i][x][0];    // �]��code rate=1/3,���codeword branch�̤j���~�t��3�A�G��3�h��
				}
			}
			else{      //�Ҽ{ i<7 �����p
				for(int x=0; x<pow(2,i+1); x++){     //x���d��̤j��2��i+1���� �Y 0<=x<2��(0~6)+1����
					k=3*i;
					for (int y=0; y<3; y++){
						if(receive[k]!=output[x][y]){    //�p�G receive���ȩM output���Ȥ��@��  weight�|�W�[�@ (output���Ȫ�ܸg��generation matrix�B��᪺��)
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

		int route[1000][128];  //���|
		for (int i=1; i<(length/3); i++){
			//�ѫe���weight�ۤ�
			if(i<7){     //�C��state�u���ߤ@�@�����|���|
				for(int x=0; x<pow(2,i); x++){
					for (int y=2*x; y<2+2*x; y++){
						weight[i][y][0]=weight[i-1][x][0]+weight[i][y][0];
						//�����ҿ�ܪ�route
						route[i][y]=x;  //��route���x�}�O���bi�ɡA�o��state�O�q�W�@��i�����ت��A�L�Ӫ�
					}
				}
			}
			else{    //����P�@state���P�ӷ���weight �ì����Ӧۭ���state
				for(int x=0; x<128; x=x+2){
					for(int y=0; y<2; y++){
						weight[i][x+y][1]=weight[i][x+y][1]+weight[i-1][x/2][0];
						weight[i][x+y][0]=weight[i][x+y][0]+weight[i-1][(x/2)+64][0];
						if(weight[i][x+y][0] > weight[i][x+y][1]){   //����U��weight
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
		//��X�̨�route�ѥXreceive code
		int back;
		int decode[1000];

		for (int i=(length/3)-1; i>0; i--){  //�q�᭱���e�s���_�ӡA��X�̨θ��|
			//printf("%d",length);
			back=route[i][back];   //�̫�@�Ӯɶ�i-1�ɡA�O�qstate0�}�l��^�h�A��Ȭ�i-2�ɪ�state
			//printf("%d",length);
			if(mem_state[back][0]!=0){  //�Y mem_state��0�A�ѽX�P�_��1
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
		//encode�p��

		int encode[3000][3];
		int memory2[7]={0};
		for(int i=0; i<(length/3); i++){     //�`�N���ɪ�encode��state2�Dstate �����ϥάۦP�� g0,g1,g2
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
		//�p��error�`��
		int error=0;     //��error�ƶq�O���~�ƶq
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
