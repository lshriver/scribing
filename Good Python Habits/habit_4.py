 # Habit 4 - Type Annotations
def upper_everything(elements: list[str]) -> list[str]:
  return [element.upper() for element in elements]

loud_list = list[str] = upper_everything(['renzo', 'jesse', 'paul', 'allison'])
