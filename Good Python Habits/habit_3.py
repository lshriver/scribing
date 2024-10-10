def is_an_adult(age: int, has_id: bool) -> bool:
  return age >= 21 and has_id

def is_bob(name: str) -> bool:
  return name.lower() =='bob'

# Habit 3: 'bundle up' above definitions into a 'mother' definition
def enter_club(name: str, age: int, has_id: bool) -> None:
  if is_bob(name):
    print('Get out of here Bob, we don\'t want no trouble.')
    return

  if is_an_adult(age, has_id):
    print('You may enter the club.')
  else:
    print('You may not enter the club.')

# Habit 2
def main() -> None:    # test cases
  enter_club('Bob', 29, has_id=True)
  enter_club('James', 29, has_id=True)
  enter_club('Sandra', 29, has_id=False)
  enter_club('Mario', 20, has_id=True)

# Habit 1
if __name__ == '__main__':
  main()
