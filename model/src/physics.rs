use cellulars_lib::constants::FloatType;
use cellulars_lib::prelude::{Pos, Rect};
use std::time::Instant;
use rand::Rng;
use cellulars_lib::lattice::Lattice;
use cellulars_lib::spin::Spin;

const GRAV_ACCEL: FloatType = -100.;
const ELASTICITY: FloatType = 0.8;

pub struct Physics {
    pub balls: Vec<Ball>,
    pub last_update: Instant
}

impl Physics {
    pub fn new(ball_radius: FloatType, n_balls: usize, rng: &mut impl Rng) -> Self {
        let mut balls = Vec::with_capacity(n_balls);
        for _ in 0..n_balls {
            balls.push(Ball::new(
                Pos::new(
                    rng.random_range(ball_radius..512. - ball_radius),
                    rng.random_range(ball_radius..424. - ball_radius),
                ),
                ball_radius,
            ))
        }
        Self {
            balls,
            last_update: Instant::now()
        }
    }

    pub fn step(&mut self, spin_latt: &Lattice<Spin>) {
        let now = Instant::now();
        let dt = (now - self.last_update).as_secs_f64();
        self.last_update = now;

        for ball in &mut self.balls {
            Self::update_ball(ball, dt, spin_latt);
        }
    }
    
    fn update_ball(ball: &mut Ball, dt: FloatType, spin_latt: &Lattice<Spin>) {
        ball.center.x += ball.vel_x * dt;
        ball.center.y += ball.vel_y * dt;
        if ball.center.y - ball.radius < 0. {
            ball.center.y = ball.radius;
            ball.vel_y = -ball.vel_y * ELASTICITY;
        }
        if ball.center.x - ball.radius < 0. {
            ball.center.x = ball.radius;
            ball.vel_x = -ball.vel_x * ELASTICITY;
        }
        if ball.center.x + ball.radius > 512. {
            ball.center.x = 512. - ball.radius;
            ball.vel_x = -ball.vel_x * ELASTICITY;
        }

        let mut closest = Pos::new(0., 0.);
        let mut closest_dist = ball.radius;
        for pos in ball.rectangle().iter_positions() {
            if pos.y >= spin_latt.rect.height() as u32 || spin_latt[pos.cast_as()] != Spin::Solid {
                continue;
            }
            let pos_f: Pos<FloatType> = pos.cast_as();
            let dx = pos_f.x - ball.center.x;
            let dy = pos_f.y - ball.center.y;
            let dist = dx.hypot(dy);
            if dist < closest_dist {
                closest_dist = dist;
                closest = pos.cast_as();
            }
        }
        if closest_dist < ball.radius {
            let dx = closest.x - ball.center.x;
            let dy = closest.y - ball.center.y;
            let ux = dx / closest_dist;
            let uy = dy / closest_dist;
            let force = (ball.radius - closest_dist) * 100. + 100.;
            // Resetting the velocities here prevents the ball from going too crazy
            ball.vel_x = -ux * force;
            ball.vel_y = -uy * force;
        }

        ball.vel_y += GRAV_ACCEL * dt;
    }
}

#[derive(Clone)]
pub struct Ball {
    pub center: Pos<FloatType>,
    pub vel_x: FloatType,
    pub vel_y: FloatType,
    pub radius: FloatType,
}

impl Ball {
    pub fn new(center: Pos<FloatType>, radius: FloatType) -> Self {
        Self {
            center,
            radius,
            vel_x: 100.,
            vel_y: 0.,
        }
    }

    pub fn rectangle(&self) -> Rect<u32> {
        Rect::new(
            Pos::new((self.center.x - self.radius).floor() as u32, (self.center.y - self.radius).floor() as u32),
            Pos::new((self.center.x + self.radius).ceil() as u32, (self.center.y + self.radius).ceil() as u32)
        )
    }
}