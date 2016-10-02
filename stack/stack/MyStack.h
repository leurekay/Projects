//#ifndef MYSTACK_H
//#define MYSTACK_H
class MyStack
{
public:
	MyStack(int size);  //分配内存，初始化栈空间
	~MyStack();          //回收栈空间内存
	bool stackEmpty(); //判断是否为空栈
	bool stackFull();
	void clearStack();   //清空栈
	int stackLength();
	//bool push(char elem); //元素入栈，栈顶上升
	//bool pop(char &elem);  //元素出栈，栈顶下降
	//void stackTraverse();    //遍历栈中所有元素

private:
	char *m_pBuffer;    //栈空间指针
	int m_iSize;				//栈容量
	int m_iTop;				//栈顶，战中元素个数
};
//#endif