
class Parent:
    
    def __init__(self, name, gender, age):
        self.name = name
        self.gender = gender
        self.age = age
        
    def __str__(self):
        
        return self.name + ' is a grown-up and has no favorite anything!'

class Child(Parent):
    def __init__(self, name, gender, age, color):
        
        super().__init__(name,gender,age)
        self.color = color
        
    def __str__(self):
        if self.gender == 'f':
            return self.name + ' is a cute little girl of ' + str(self.age) + ' years!'+ '\n' + 'Her favorite color is ' + self.color
            
        elif self.gender == 'm':
            return self.name + ' is a strong lad of '+ str(self.age) + ' years!' + '\n' + 'His favorite color is ' + self.color
            
        

Steve = Parent('Steve','m', 32)
Eva = Child('Eva','f', 5, 'pink')
John = Child('John','m',7,'blue')

print(Eva)
print(John)
print(Steve)