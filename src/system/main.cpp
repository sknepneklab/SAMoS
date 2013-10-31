#include "system.hpp"
#include "../messenger/messenger.hpp"

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
  
  string name = "meho";

  MessengerPtr msg( new Messenger(name));

  SystemPtr sys(new System("data.dat",msg));

  cout << sys->get_particle(1);

  Particle& p = sys->get_particle(1);

  p.x += 1.0;

  cout << sys->get_particle(1);

    
  return 0;
}
