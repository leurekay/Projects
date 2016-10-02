//#ifndef MYSTACK_H
//#define MYSTACK_H
class MyStack
{
public:
	MyStack(int size);  //�����ڴ棬��ʼ��ջ�ռ�
	~MyStack();          //����ջ�ռ��ڴ�
	bool stackEmpty(); //�ж��Ƿ�Ϊ��ջ
	bool stackFull();
	void clearStack();   //���ջ
	int stackLength();
	//bool push(char elem); //Ԫ����ջ��ջ������
	//bool pop(char &elem);  //Ԫ�س�ջ��ջ���½�
	//void stackTraverse();    //����ջ������Ԫ��

private:
	char *m_pBuffer;    //ջ�ռ�ָ��
	int m_iSize;				//ջ����
	int m_iTop;				//ջ����ս��Ԫ�ظ���
};
//#endif