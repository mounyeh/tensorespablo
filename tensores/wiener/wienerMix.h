
struct vtkFileHeader {
	char fileName[100];
	unsigned dimensions[3];
	float spacing[3];
	float origin[3];
	unsigned size;
	float *data;
};

void vtkFileHeaderReader(struct vtkFileHeader *);

void vtkFileDataReader(struct vtkFileHeader *);

void vtkFileHeaderWriter(struct vtkFileHeader);

void vtkFileDataWriter(struct vtkFileHeader);

void filterData(struct vtkFileHeader,struct vtkFileHeader *,float,float, unsigned, unsigned, unsigned);

void swap_4(void *);
